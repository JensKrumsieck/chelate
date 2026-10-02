use std::{fmt, io};
use thiserror::Error;

use crate::types::InvalidElementError;

/// What went wrong while reading a file.
#[derive(Debug, Error)]
pub enum ErrorKind {
    #[error(transparent)]
    InvalidElement(#[from] InvalidElementError),
    #[error(transparent)]
    Io(#[from] io::Error),
    /// Input that does not match the expected file format.
    #[error("{0}")]
    Parse(String),
}

/// Error while reading a file, with the 1-based line it belongs to if known.
#[derive(Debug)]
pub struct FileError {
    pub kind: ErrorKind,
    pub line_number: Option<usize>,
}

impl FileError {
    /// Creates a format error that does not belong to a specific line.
    /// # Examples
    /// ```
    /// use chelate::error::FileError;
    ///
    /// let error = FileError::parse("Invalid counts line");
    ///
    /// assert_eq!(error.line_number, None);
    /// assert_eq!(error.to_string(), "Invalid counts line");
    /// ```
    pub fn parse(message: impl Into<String>) -> Self {
        ErrorKind::Parse(message.into()).into()
    }

    /// Sets the 1-based number of the line the error belongs to.
    /// Line parsers usually do not know their line number, so it can be added where the lines are iterated.
    /// # Examples
    /// ```
    /// use chelate::error::FileError;
    ///
    /// let error = FileError::parse("Invalid counts line").with_line(4);
    ///
    /// assert_eq!(error.line_number, Some(4));
    /// assert_eq!(error.to_string(), "line 4: Invalid counts line");
    /// ```
    pub fn with_line(mut self, line_number: usize) -> Self {
        self.line_number = Some(line_number);
        self
    }
}

impl fmt::Display for FileError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.line_number {
            Some(line_number) => write!(f, "line {}: {}", line_number, self.kind),
            None => self.kind.fmt(f),
        }
    }
}

impl std::error::Error for FileError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        self.kind.source()
    }
}

impl From<ErrorKind> for FileError {
    fn from(kind: ErrorKind) -> Self {
        FileError {
            kind,
            line_number: None,
        }
    }
}

impl From<io::Error> for FileError {
    fn from(error: io::Error) -> Self {
        ErrorKind::from(error).into()
    }
}

impl From<InvalidElementError> for FileError {
    fn from(error: InvalidElementError) -> Self {
        ErrorKind::from(error).into()
    }
}
