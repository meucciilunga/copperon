use crate::search;
use crate::analysis;
use std::fs;
use std::path;
use std::sync::Mutex;

trait Log {
    fn log_summary_entry(&self, shared_log_file: &Mutex<fs::File>);
}

