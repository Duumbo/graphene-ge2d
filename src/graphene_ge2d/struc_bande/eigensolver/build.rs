use pkg_config;

fn main() {
    pkg_config::Config::new().probe("lapack").unwrap();
    println!("cargo:rerun-if-changed=build.rs");
}
