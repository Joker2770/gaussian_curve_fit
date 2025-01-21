use thiserror::Error;

#[derive(Debug, Error)]
pub enum CustomError {
    #[error("Illegal sigma(sigma: {0})")] 
    IllegalSigma(f32),
    #[error("Illegal sigma(sigma_x: {0}, sigma_y: {1})")]
    IllegalSigmaXY(f32, f32),
    #[error("Illegal result from SVD-solve")]
    IllegalResult,
    #[error("SVD solve failed")]
    SvdFailed,
    #[error("Matrix size not match")]
    MatrixSizeNotMatch,
}
