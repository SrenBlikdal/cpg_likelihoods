use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

static NEXT_ID: AtomicU64 = AtomicU64::new(0);

pub struct Scratch(PathBuf);

impl Scratch {
    pub fn new() -> Self {
        let id = NEXT_ID.fetch_add(1, Ordering::Relaxed);
        let root = std::env::temp_dir().join(format!(
            "cpg_likelihoods_test_{}_{}_{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos(),
            id
        ));
        std::fs::create_dir(&root).unwrap();
        Self(root)
    }

    pub fn path(&self, name: &str) -> PathBuf {
        self.0.join(name)
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

pub fn write_input(path: &Path) {
    std::fs::write(
        path,
        concat!(
            "chr1\t0\t1\tC\t30\t+\t0\t1\t255,0,0\t30\t0\t0\t30\t0\t0\t0\t0\t2\n",
            "chr1\t1\t2\tC\t15\t+\t1\t2\t255,0,0\t15\t0\t0\t15\t0\t0\t0\t15\t1\n",
            "chr1\t2\t3\tC\t1\t+\t2\t3\t255,0,0\t1\t0\t0\t1\t0\t0\t0\t28\t1\n",
        ),
    )
    .unwrap();
}
