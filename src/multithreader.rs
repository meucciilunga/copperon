#![allow(dead_code)]
use std::sync::mpsc::Receiver;
use std::sync::{Arc, Mutex, mpsc};
use std::thread;

enum DataSignal {
    Data(ThreadPoolTask),
    Exit,
}

enum WorkerSignal {
    Idle,
    Busy,
    Closed,
}

pub struct ThreadPoolTask {
    task_id: usize,
    subroutine: Box<dyn FnOnce() + Send>,
}

impl ThreadPoolTask {
    fn new(task_id: usize, subroutine: Box<dyn FnOnce() + Send + 'static>) -> ThreadPoolTask {
        ThreadPoolTask { 
            task_id, 
            subroutine,
        }
    }
}

pub struct Manager {
    sequential_lock_tx: mpsc::Sender<DataSignal>,     // Universal sequential transmitter
    status_update_rx: mpsc::Receiver<WorkerSignal>,
    num_workers: usize,
}

impl Manager {

    pub fn new(num_threads: usize) -> Manager {
        
        // Limitations on number of threads
        assert!(num_threads > 0);
        assert!(num_threads <= 200);

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

            // Activate worker and make their thread available to manager
            let new_worker = Worker::new(thread_id, tx, rx);
            new_worker.activate();
        }

        // Construct manager
        Manager {
            status_update_rx: status_rx,
            sequential_lock_tx: seq_task_tx,
            num_workers: num_threads,
        }
    }

    fn close_threadpool(&self) {
        for _ in 0..self.num_workers {
            self.sequential_lock_tx.send(DataSignal::Exit).expect("ERROR: could not send exit signal to threadpool!");
        }
        
        let mut closed_count = 0_usize;
        loop {

            // Exit loop when all created worker threads have been CONFIRMED as closed
            if closed_count == self.num_workers {
                break
            }

            // Pull backlogged responses from update rx
            let receipt_err = "ERROR: while closing threadpool, could not receive 'Closed' signal from pool.";
            let receipt = self.status_update_rx.recv().expect(receipt_err);

            match receipt {
                WorkerSignal::Closed => closed_count += 1,
                _ => continue,
            }
        }
    }

    pub fn execute_queue(&self, mut task_list: Vec<ThreadPoolTask>) {

        // While queue still has tasks available
        while task_list.len() > 0 {
            let receipt_err = "ERROR: while distributing tasks to threadpool, could not receive status update from pool.";
            let receipt = self.status_update_rx.recv().expect(receipt_err);

            if let WorkerSignal::Idle = receipt {
                let new_task = task_list.pop();

                if let Some(task) = new_task {
                    let send_err = "ERROR: management thread could not send task to threadpool!";
                    let task = DataSignal::Data(task);
                    self.sequential_lock_tx.send(task).expect(send_err);
                }
            }
        }

        // Send signal to close threadpool after all tasks have been distributed
        self.close_threadpool();
    }
}

struct Worker {
    id: usize,
    data_rx: Option<Arc<Mutex<Receiver<DataSignal>>>>,
    status_tx: Option<mpsc::Sender<WorkerSignal>>,
}

impl Worker {
    fn new(id: usize, 
           thread_status_tx: mpsc::Sender<WorkerSignal>, 
           thread_data_rx: Arc<Mutex<Receiver<DataSignal>>>) -> Worker
        {
            Worker { 
                id,
                data_rx: Some(thread_data_rx), 
                status_tx: Some(thread_status_tx), 
            }
    }

    fn activate(mut self) {
        
        // Allow local function variables to obtain ownership of a copies of respective channel ends
        let thread_tx = self.status_tx.expect("ERROR: worker did not have channel Status Tx available for spawning a new thread.");
        let thread_rx = self.data_rx.expect("ERROR: worker did not have channel DataSignal Rx available for spawning a new thread.");
        let thread_id = self.id;

        // Mark worker tx/rx fields as empty
        self.status_tx = None;
        self.data_rx = None;

        // List of error messages
        let lock_err = format!("ERROR: could not acquire lock for thread_{}!", thread_id);
        let idle_err = format!("ERROR: thread_{} could not deliver 'Idle' signal flag to management thread!", thread_id);
        let receipt_err = format!("ERROR: thread_{} could not receive task from the management thread! Channel hung up.", thread_id);
        let busy_err = format!("ERROR: thread_{} could not deliver 'Busy' signal flag to management thread!", thread_id);
        let closed_err = format!("ERROR: thread_{} could not deliver 'Closed' signal flag to management thread!", thread_id);
        let spawn_err = format!("ERROR: thread_{} could not be spawned!", thread_id);

        // Build thread loop logic + import the communication
        // channel elements taken from original worker struct
        let open_thread = move || {
            
            // Force thread to idle until a task from the queue is available; run thread
            // until an 'Exit' signal is received.
            loop {
                // Change worker status to available and wait to acquire a DataSignal from manager
                thread_tx.send(WorkerSignal::Idle).expect(&idle_err);

                // Attempt to acquire lock on queue data; once acquired, receive data
                let task_receipt = thread_rx.lock().expect(&lock_err).recv().expect(&receipt_err);

                //  Exit loop or pull wrapped closure for execution
                let task = match task_receipt {
                    DataSignal::Data(s) => s,
                    DataSignal::Exit => break,
                };

                // Change worker status and run submitted task
                thread_tx.send(WorkerSignal::Busy).expect(&busy_err);
                (task.subroutine)();
            }

            // Change worker status to closed and allow thread to shutdown
            thread_tx.send(WorkerSignal::Closed).expect(&closed_err);
            println!("thread_{} has been closed.", thread_id);
        };

        // Formally instantiate the new thread w/ cutstom name
        let worker_thread = thread::Builder::new().name(thread_id.to_string());
        worker_thread.spawn(open_thread).expect(&spawn_err);
    }
}

#[cfg(test)]
mod test {
    use std::time::{SystemTime, Duration};
    use super::*;

    #[test]
    #[should_panic]
    fn open_close_pool_1() {
        let test_threadpool = Manager::new(0);
        test_threadpool.close_threadpool();
    }

    #[test]
    #[should_panic]
    fn open_close_pool_2() {
        let test_threadpool = Manager::new(314);
        test_threadpool.close_threadpool();
    }

    #[test]
    fn open_close_pool_3() {
        let test_threadpool = Manager::new(24);
        test_threadpool.close_threadpool();
    }

    #[test]
    fn task_queue_example_1() {
        let num_threads = 24;
        let test_task_size = num_threads*10;
        let new_pool = Manager::new(num_threads);
        let mut task_list: Vec<ThreadPoolTask> = Vec::with_capacity(test_task_size);

        for id in 0..test_task_size {
            let x = Duration::from_millis(100);
            let id_copy = id;
            let new_subroutine = move || {
                thread::sleep(x);
            };
            let new_subroutine = Box::new(new_subroutine);
            let new_task = ThreadPoolTask::new(id_copy, new_subroutine);
            task_list.push(new_task);
        }

        let now = SystemTime::now();
        new_pool.execute_queue(task_list);
        let elapsed_duration = now.elapsed().unwrap();
        println!("Number of threads: {}\nNumber of tasks: {}\nTime to Process: {}s\n", num_threads, test_task_size, elapsed_duration.as_secs_f32());
    }

    #[test]
    fn task_queue_example_2() {
        let num_threads = 12;
        let test_task_size = num_threads*10*2;
        let new_pool = Manager::new(num_threads);
        let mut task_list: Vec<ThreadPoolTask> = Vec::with_capacity(test_task_size);

        for id in 0..test_task_size {
            let x = Duration::from_millis(100);
            let id_copy = id;
            let new_subroutine = move || {
                thread::sleep(x);
            };
            let new_subroutine = Box::new(new_subroutine);
            let new_task = ThreadPoolTask::new(id_copy, new_subroutine);
            task_list.push(new_task);
        }

        let now = SystemTime::now();
        new_pool.execute_queue(task_list);
        let elapsed_duration = now.elapsed().unwrap();
        println!("Number of threads: {}\nNumber of tasks: {}\nTime to Process: {}s\n", num_threads, test_task_size, elapsed_duration.as_secs_f32());
    }

    #[test]
    fn task_queue_example_3() {
        let num_threads = 3;
        let test_task_size = num_threads*10*2*4;
        let new_pool = Manager::new(num_threads);
        let mut task_list: Vec<ThreadPoolTask> = Vec::with_capacity(test_task_size);

        for id in 0..test_task_size {
            let x = Duration::from_millis(100);
            let id_copy = id;
            let new_subroutine = move || {
                thread::sleep(x);
            };
            let new_subroutine = Box::new(new_subroutine);
            let new_task = ThreadPoolTask::new(id_copy, new_subroutine);
            task_list.push(new_task);
        }

        let now = SystemTime::now();
        new_pool.execute_queue(task_list);
        let elapsed_duration = now.elapsed().unwrap();
        println!("Number of threads: {}\nNumber of tasks: {}\nTime to Process: {}s\n", num_threads, test_task_size, elapsed_duration.as_secs_f64());
    }
}