use std::sync::atomic::{AtomicUsize, Ordering};

pub mod xyz;
pub mod pdb;

pub fn normalize_symbol(symbol: &str) -> String {
    let normalized_symbol = if let Some(first_char) = symbol.chars().next() {
        first_char.to_uppercase().collect::<String>() + &symbol[1..].to_lowercase()
    } else {
        String::new()
    };
    normalized_symbol
}

pub fn get_next_id() -> usize {
    static COUNTER: AtomicUsize = AtomicUsize::new(0);
    COUNTER.fetch_add(1, Ordering::Relaxed)
}