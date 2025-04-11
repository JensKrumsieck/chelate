use std::sync::atomic::{AtomicUsize, Ordering};

pub mod cif;
pub mod mol;
pub mod mol2;
pub mod pdb;
pub mod xyz;

pub fn normalize_symbol(symbol: &str) -> String {
    let normalized_symbol = if let Some(first_char) = symbol.chars().next() {
        first_char.to_uppercase().collect::<String>() + &symbol[1..].to_lowercase()
    } else {
        String::new()
    };
    normalized_symbol
}

static COUNTER: AtomicUsize = AtomicUsize::new(0);
pub fn get_next_id() -> usize {
    COUNTER.fetch_add(1, Ordering::Relaxed)
}

pub(crate) fn reset_counter() {
    COUNTER.store(0, Ordering::Relaxed);
}
