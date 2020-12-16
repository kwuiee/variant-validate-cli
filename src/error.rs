pub fn err(v: &str) -> std::io::Error {
    std::io::Error::new(std::io::ErrorKind::InvalidData, v)
}

pub fn opterr() -> std::io::Error {
    std::io::Error::new(std::io::ErrorKind::InvalidData, "Option error.")
}
