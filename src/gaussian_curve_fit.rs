/*
A no_std and no heap memory library for gaussian curve coefficents calculation.
Copyright (C) 2024  joker2770

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#![allow(unused_imports)]

pub mod gaussian_curve {
    use core::f32;

    use nalgebra::{ComplexField, SMatrix};

    const POINT_NUM_2D: usize = 8;
    const POINT_NUM_3D: usize = 16;

    const MATRIX_COLUMN_2D: usize = 3;
    const MATRIX_COLUMN_3D: usize = 5;

    #[derive(Debug, Default)]
    pub struct GaussianCoefficents2D {
        alpha: f32,
        mu: f32,
        sigma: f32,
    }

    impl GaussianCoefficents2D {
        /// $$f(x) & = {\alpha}{e^{-{\frac {(x - \mu)^2}{2\sigma^2}}}}$$
        pub fn value(&self, x: f32) -> f32 {
            if self.sigma.abs() > 0.0f32 {
                self.alpha
                    * (-(x - self.mu).powf(2.0f32) / (2.0f32 * self.sigma * self.sigma)).exp()
            } else {
                0.0f32
            }
        }

        pub fn coefficents(&self) -> (f32, f32, f32) {
            (self.alpha, self.mu, self.sigma)
        }

        pub fn peak(&self) -> (f32, f32) {
            (self.mu, self.alpha)
        }

        pub fn get_matrix_data_from_8_points(
            x: &[f32; POINT_NUM_2D],
            y: &[f32; POINT_NUM_2D],
        ) -> ([f32; POINT_NUM_2D * MATRIX_COLUMN_2D], [f32; POINT_NUM_2D]) {
            let mut tmp_x = [0.0f32; POINT_NUM_2D * MATRIX_COLUMN_2D];
            let mut tmp_y = [0.0f32; POINT_NUM_2D];
            for i in 0..POINT_NUM_2D {
                tmp_x[(i * MATRIX_COLUMN_2D + 0) as usize] = 1.0f32;
                tmp_x[(i * MATRIX_COLUMN_2D + 1) as usize] = x[i] as f32;
                tmp_x[(i * MATRIX_COLUMN_2D + 2) as usize] = (x[i] as f32).powf(2.0f32);
                tmp_y[i] = (y[i] as f32).ln();
            }
            (tmp_x, tmp_y)
        }

        /// Solve ax = b.
        /// Any singular value smaller than `eps` is assumed to be zero.
        pub fn get_coefficents_from_8_matrix_data(
            &mut self,
            a: &[f32; POINT_NUM_2D * MATRIX_COLUMN_2D],
            b: &[f32; POINT_NUM_2D],
            eps: f32,
        ) -> Self {
            if a.len() == MATRIX_COLUMN_2D * b.len() {
                type MatrixXx1f32 = SMatrix<f32, POINT_NUM_2D, 1>;
                type MatrixXx3f32 = SMatrix<f32, POINT_NUM_2D, MATRIX_COLUMN_2D>;
                let mb = MatrixXx1f32::from_row_slice(b);
                let ma = MatrixXx3f32::from_row_slice(a);

                let decomp = ma.svd(true, true);

                let x = decomp.solve(&mb, eps);
                if let Ok(r) = x {
                    if r[2] < 0.0f32 {
                        self.alpha = (r[0] - r[1].powf(2.0f32) / (4.0f32 * r[2])).exp();
                        self.mu = -r[1] / (2.0f32 * r[2]);
                        self.sigma = 1.0f32 / (-2.0f32 * r[2]).sqrt();
                    } else {
                        self.alpha = 0.0f32;
                        self.mu = 0.0f32;
                        self.sigma = 0.0f32;
                    }
                } else {
                    self.alpha = 0.0f32;
                    self.mu = 0.0f32;
                    self.sigma = 0.0f32;
                }
            }

            Self {
                alpha: self.alpha,
                mu: self.mu,
                sigma: self.sigma,
            }
        }
    }

    #[derive(Debug, Default)]
    pub struct GaussianCoefficents3D {
        alpha: f32,
        mu_x: f32,
        mu_y: f32,
        sigma_x: f32,
        sigma_y: f32,
    }

    impl GaussianCoefficents3D {
        /// $$f(x, y) & = {\alpha}{e^{{-{\frac {(x - {\mu}_x)^2}{2{\sigma}_x^2}}} + {-{\frac {(y - {\mu}_y)^2}{2{\sigma}_y^2}}}}}$$
        pub fn value(&self, x: f32, y: f32) -> f32 {
            if self.sigma_x.abs() > 0.0f32 && self.sigma_y.abs() > 0.0f32 {
                self.alpha
                    * ((-(x - self.mu_x).powf(2.0f32) / (2.0f32 * self.sigma_x.powf(2.0f32)))
                        + (-(y - self.mu_y).powf(2.0f32) / (2.0f32 * self.sigma_y.powf(2.0f32))))
                    .exp()
            } else {
                0.0f32
            }
        }

        pub fn peak(&self) -> (f32, f32, f32) {
            (self.mu_x, self.mu_y, self.alpha)
        }

        pub fn coefficents(&self) -> (f32, f32, f32, f32, f32) {
            (self.alpha, self.mu_x, self.mu_y, self.sigma_x, self.sigma_y)
        }

        pub fn get_matrix_data_from_16_points(
            x: &[f32; POINT_NUM_3D],
            y: &[f32; POINT_NUM_3D],
            z: &[f32; POINT_NUM_3D],
        ) -> ([f32; POINT_NUM_3D * MATRIX_COLUMN_3D], [f32; POINT_NUM_3D]) {
            let mut tmp_x_y = [0.0f32; POINT_NUM_3D * MATRIX_COLUMN_3D];
            let mut tmp_z = [0.0f32; POINT_NUM_3D];
            for i in 0..POINT_NUM_3D {
                tmp_x_y[(i * MATRIX_COLUMN_3D + 0) as usize] = 1.0f32;
                tmp_x_y[(i * MATRIX_COLUMN_3D + 1) as usize] = x[i] as f32;
                tmp_x_y[(i * MATRIX_COLUMN_3D + 2) as usize] = y[i] as f32;
                tmp_x_y[(i * MATRIX_COLUMN_3D + 3) as usize] = (x[i] as f32).powf(2.0f32);
                tmp_x_y[(i * MATRIX_COLUMN_3D + 4) as usize] = (y[i] as f32).powf(2.0f32);
                tmp_z[i] = (z[i] as f32).ln();
            }
            (tmp_x_y, tmp_z)
        }

        /// Solve ax = b.
        /// Any singular value smaller than `eps` is assumed to be zero.
        pub fn get_coefficents_from_16_matrix_data(
            &mut self,
            a: &[f32; POINT_NUM_3D * MATRIX_COLUMN_3D],
            b: &[f32; POINT_NUM_3D],
            eps: f32,
        ) -> Self {
            if a.len() == MATRIX_COLUMN_3D * b.len() {
                type MatrixXx1f32 = SMatrix<f32, POINT_NUM_3D, 1>;
                type MatrixXx5f32 = SMatrix<f32, POINT_NUM_3D, MATRIX_COLUMN_3D>;
                let mb = MatrixXx1f32::from_row_slice(b);
                let ma = MatrixXx5f32::from_row_slice(a);

                let decomp = ma.svd(true, true);

                let x = decomp.solve(&mb, eps);
                if let Ok(r) = x {
                    if r[3] < 0.0f32 && r[4] < 0.0f32 {
                        self.alpha = (r[0]
                            - r[1].powf(2.0f32) / (4.0f32 * r[3])
                            - r[2].powf(2.0f32) / (4.0f32 * r[4]))
                            .exp();
                        self.mu_x = -r[1] / (2.0f32 * r[3]);
                        self.mu_y = -r[2] / (2.0f32 * r[4]);
                        self.sigma_x = 1.0f32 / (-2.0f32 * r[3]).sqrt();
                        self.sigma_y = 1.0f32 / (-2.0f32 * r[4]).sqrt();
                    } else {
                        self.alpha = 0.0f32;
                        self.mu_x = 0.0f32;
                        self.mu_y = 0.0f32;
                        self.sigma_x = 0.0f32;
                        self.sigma_y = 0.0f32;
                    }
                } else {
                    self.alpha = 0.0f32;
                    self.mu_x = 0.0f32;
                    self.mu_y = 0.0f32;
                    self.sigma_x = 0.0f32;
                    self.sigma_y = 0.0f32;
                }
            }

            Self {
                alpha: self.alpha,
                mu_x: self.mu_x,
                mu_y: self.mu_y,
                sigma_x: self.sigma_x,
                sigma_y: self.sigma_y,
            }
        }
    }
}
