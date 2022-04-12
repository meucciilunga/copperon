use std::sync::{mpsc, Arc, Mutex};
use std::thread;
use std::fs;
use std::path::PathBuf;
use crate::genome_search::*;


pub struct GenomeTask {
    task_id: String,
    subroutine: Box<dyn FnOnce() + Send + 'static>,
}

// Build a closure that can process entire genome
impl GenomeTask {
    // Builds a closure that can process an entire genome start to finish
    pub fn new(name: String,
           genome_path: String,
           annotation_path: String,
           search_rad: usize) -> GenomeTask
        
        {
            let id = genome_path.clone();

            // Closure to run
            let processing_task = move || {

                // Consensus sequence to search for
                let operator_seq = "RNYKACANNYGTMRNY";

                // Build genome object
                let mut genome = CompleteGenome::new_annotated(name,
                                                               &genome_path[..], 
                                                               &annotation_path[..]);

                // Analyze genome object for consensus sequence matches
                genome.analyze(operator_seq, search_rad);

                // Print results to file (by first defining an output dir path)
                let mut output_dir_path = PathBuf::new();
                output_dir_path.push("/home/katanga/Desktop/test_downloads/test_blast");

                genome.export_results(output_dir_path);
            };
            
            GenomeTask { 
                task_id: id, 
                subroutine: Box::new(processing_task),
            }
            
    }

    // Build new GenomeTask when given the location of a genome subdir
    pub fn task_from_subdir(subdir: String, search_rad: usize) -> GenomeTask {

        // Pull name from subdirectory path
        let name = std::path::Path::new(&subdir);
        let name = name.file_name().unwrap().to_str().unwrap().to_string();

        // Define suffixes that each needed file uses; each file is named
        // after the name of the subdirectory, plus one of these prefixes
        let mut genome_file = "_genomic.fna".to_string();
        let mut feature_table = "_feature_table.txt".to_string();
        let mut assembly_report = "_assembly_report.txt".to_string();

        // Read in all files from given genome subdirectory
        let paths = fs::read_dir(subdir).unwrap();

        // Utilize file based on the suffix it contains
        for path in paths {
            let curr_path = path.unwrap().path().to_str().unwrap().to_string();
            
            if curr_path.contains(&genome_file) {
                genome_file = curr_path;
            }
            else if curr_path.contains(&feature_table) {
                feature_table = curr_path;
            }
            else if curr_path.contains(&assembly_report) {
                assembly_report = curr_path;
            }
        }

        GenomeTask::new(name, genome_file.to_string(), feature_table.to_string(), search_rad)
    }
}

enum DataSignal {
    Data(GenomeTask),
    Exit,
}

enum WorkerSignal {
    Idle(usize),
    Busy(usize),
    Closed(usize),
}

// Pool Architecture
pub struct Manager {
    pub workers: Vec<Worker>,
    pub num_active_workers: usize,
    sequential_lock_tx: mpsc::Sender<DataSignal>,     // Universal sequential transmitter
    status_update_rx: mpsc::Receiver<WorkerSignal>,
}

impl Manager {

    // Build new manager
    pub fn new(num_threads: usize) -> Manager {

        // Can't have a pool with zero threads
        assert!(num_threads > 0);
        assert!(num_threads <= 200);

        // Pre-allocate space for a vector of workers
        let mut workers: Vec<Worker> = Vec::with_capacity(num_threads);

        // Build global communication infrastructure
        let (seq_task_tx, seq_task_rx): (mpsc::Sender<DataSignal>, mpsc::Receiver<DataSignal>) = mpsc::channel();
        let (status_tx, status_rx): (mpsc::Sender<WorkerSignal>, mpsc::Receiver<WorkerSignal>) = mpsc::channel();
        let lockable_rx = Mutex::new(seq_task_rx);
        let seq_lock_rx = Arc::new(lockable_rx);
        // NOTE: original 'Arc' should be dropped at the end of call to new(),
        // leaving--at the end of this function call--only as many references 
        // to lockable_rx as there are instantiated worker threads
        
        // Build same number of workers as threads_requested
        for thread_id in 0..num_threads {

            // Give workers copies of their end of the communicaiton infrastructure
            let rx = Arc::clone(&seq_lock_rx);
            let tx = status_tx.clone();

            // Add worker to the manager's  list
            let new_worker = Worker::new(thread_id, Some(rx), Some(tx));
            workers.push(new_worker);
        }
        
        Manager {
            workers: workers,
            num_active_workers: 0,
            sequential_lock_tx: seq_task_tx,
            status_update_rx: status_rx,
        }
    }

    // Activate worker threads
    fn activate_workers(&mut self) {
        for worker in &mut self.workers {
            self.num_active_workers += 1;
            worker.activate();
        }
    }

    // Have threadpool manager close an active thread
    fn kill_worker(&mut self) {
        let task = DataSignal::Exit;
        self.send_to_pool(task);
        self.num_active_workers -= 1;
    }

    // Have manager send a DataSignal out to the task channel
    fn send_to_pool(&mut self, task: DataSignal) {
        let err_msg = "ERROR: Manager could not reach worker pool for task submission!";
        self.sequential_lock_tx.send(task).expect(err_msg);
    }

    // Have threadpool manager distribute a single task to pool for processing
    fn execute(&mut self, job: GenomeTask) {
        let task = DataSignal::Data(job);
        self.send_to_pool(task);
    }

    // Have threadpool manager distribute a list of tasks
    pub fn execute_queue(&mut self, mut task_list: Vec<GenomeTask>) {

        self.activate_workers();
        let mut closed_count = 0_usize;
        
        // Only hand out jobs as workers become available
        while task_list.len() != 0 {

            // Check for an availability signal; over the lifetime of the program, you'll 
            // always have more availability flags than tasks, preventing deadlock
            match self.status_update_rx.recv() { // Pull a message from backlog or wait until one arrives
                Ok(WorkerSignal::Idle(_)) => {
                    let task = task_list.pop().unwrap();
                    self.execute(task); 
                },
                Ok(WorkerSignal::Busy(_)) => {
                    continue;
                },
                _ => { panic!("ERROR: job distribution subroutine received 'Closed' flag.") }
            }
        }

        // Start closing (idle) threads after last job is distributed
        while closed_count < self.workers.len() {
            
            match self.status_update_rx.recv() { // Pull a message from backlog or wait until one arrives

                // Any remaining idle_flags indicate a thread can be closed since no more tasks remain
                Ok(WorkerSignal::Idle(_)) => {
                    self.kill_worker();
                },

                // Don't count a thread as closed until its "Closed" flag is received
                Ok(WorkerSignal::Closed(_)) => {
                    closed_count += 1;
                },
                // this ensures the manager doesn't kill the program before the last
                // thread has truly finished its work

                _ => (),  // Ignore busy flags
            }
        }
    }
}

pub struct Worker {
    id: usize,
    sequential_lock_rx: Option<Arc<Mutex<mpsc::Receiver<DataSignal>>>>,  // Universal sequential receiver
    status_update_tx: Option<mpsc::Sender<WorkerSignal>>,
}

impl Worker {

    // Build new worker
    fn new(id: usize, 
           sequential_lock_rx: Option<Arc<Mutex<mpsc::Receiver<DataSignal>>>>,
           status_update_tx: Option<mpsc::Sender<WorkerSignal>>) -> Worker {
        
        Worker {
            id,
            sequential_lock_rx,
            status_update_tx,
        }
    }

    // Activate new worker's thread
    fn activate(&mut self) {

        // Gain ownership over a copy of each of the worker's channel ends
        let worker_tx = self.status_update_tx.clone().expect("ERROR: worker has no status tx available for use in spawned thread.");
        let worker_rx = self.sequential_lock_rx.clone().expect("ERROR: worker has no sequential rx available for use in spawned thread.");

        // Clone the channel ends using respective methods (+ id)
        let thread_tx = worker_tx.clone();
        let thread_rx = Arc::clone(&worker_rx);
        let thread_id = self.id; // So manager can know which worker is signalling

        // Mark worker tx/rx fields as empty
        self.status_update_tx = None;
        self.sequential_lock_rx = None;

        // List of error messages
        let msg = format!("ERROR: could not acquire lock for thread_{}!", thread_id);
        let idle_err = format!("thread_{} could not deliver 'Idle' signal flag to management thread!", thread_id);
        let busy_err = format!("thread_{} could not deliver 'Busy' signal flag to management thread!", thread_id);
        let closed_err = format!("thread_{} could not deliver 'Closed' signal flag to management thread!", thread_id);
        let spawn_err = format!("thread_{} could not be spawned!", thread_id);
        
        // Build thread loop logic + import the copies of the communication 
        // channel elements from original worker struct
        let open_thread = move || {
            
            // Force thread to idle until a task from the queue is available; run thread
            // until an 'Exit' signal is received.
            loop {
                // Change worker status to available and wait to acquire a task
                thread_tx.send(WorkerSignal::Idle(thread_id)).expect(&idle_err);

                // Attempt to acquire lock on queue data
                let task_receipt = thread_rx.lock().expect(&msg).recv().unwrap();

                //  Exit loop or pull wrapped closure for execution
                let task = match task_receipt {
                    DataSignal::Data(s) => { s },
                    DataSignal::Exit => { break }
                };

                // Change worker status and run submitted task
                thread_tx.send(WorkerSignal::Busy(thread_id)).expect(&busy_err);
                println!("Worker #{}: executing task '{}'", thread_id, task.task_id);
                (task.subroutine)();
            }

            // Change worker status to closed and allow thread to shutdown
            thread_tx.send(WorkerSignal::Closed(thread_id)).expect(&closed_err);
            println!("thread_{} has been closed.", thread_id);
        };

        // Formally instantiate the new thread w/ cutstom name
        let worker_thread = thread::Builder::new().name(thread_id.to_string());
        worker_thread.spawn(open_thread).expect(&spawn_err);

        // DEBUG OUTPUT
        println!("Worker #{} has been activated!", self.id);
    }

}