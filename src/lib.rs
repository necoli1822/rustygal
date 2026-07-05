// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2007-2016 Doug Hyatt, Univ. of Tennessee / UT-Battelle (original Prodigal C)
// Copyright (C) 2026 Sunju Kim (Rust reimplementation)

pub mod bitmap;
pub mod dprog;
pub mod fptr;
pub mod gene;
pub mod metagenomic;
pub mod node;
pub mod sequence;
pub mod training;
pub mod translation;

pub mod api;
pub mod meta_api;

#[cfg(feature = "python")]
pub mod python;

pub use bitmap::*;
pub use dprog::*;
pub use gene::*;
pub use metagenomic::*;
pub use node::*;
pub use sequence::*;
pub use training::*;
pub use translation::*;

pub const VERSION: &str = "0.2.0";
pub const DATE: &str = "June, 2026";
