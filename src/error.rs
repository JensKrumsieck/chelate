use std::{error::Error, fmt, io};
use thiserror::Error;

#[derive(Debug, Error)]
pub enum FileError {
    #[error(transparent)]
    Io(#[from] io::Error),
    #[error(transparent)]
    Parse(#[from] ParseError),
}

/// Error for input that does not match the expected file format.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ParseError {
    pub message: String,
    pub line_number: Option<usize>,
}

impl ParseError {
    /// Creates an error that does not belong to a specific line.
    /// # Examples
    /// ```
    /// use chelate::error::ParseError;
    ///
    /// let error = ParseError::new("Invalid counts line");
    ///
    /// assert_eq!(error.line_number, None);
    /// assert_eq!(error.to_string(), "Invalid counts line");
    /// ```
    pub fn new(message: impl Into<String>) -> Self {
        ParseError {
            message: message.into(),
            line_number: None,
        }
    }

    /// Sets the 1-based number of the line the error belongs to.
    /// Line parsers usually do not know their line number, so it can be added where the lines are iterated.
    /// # Examples
    /// ```
    /// use chelate::error::ParseError;
    ///
    /// let error = ParseError::new("Invalid counts line").with_line(4);
    ///
    /// assert_eq!(error.line_number, Some(4));
    /// assert_eq!(error.to_string(), "line 4: Invalid counts line");
    /// ```
    pub fn with_line(mut self, line_number: usize) -> Self {
        self.line_number = Some(line_number);
        self
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.line_number {
            Some(line_number) => write!(f, "line {}: {}", line_number, self.message),
            None => f.write_str(&self.message),
        }
    }
}

impl Error for ParseError {}
