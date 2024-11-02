/*
A no_std library for gaussian curve coefficents calculation.
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

    const POINT_CNT: usize = 8;

    #[derive(Debug, Default)]
    pub struct GaussianCoefficents {
        alpha: f32,
        mu: f32,
        sigma: f32,
    }

    impl GaussianCoefficents {
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
            x: &[f32; POINT_CNT],
            y: &[f32; POINT_CNT],
        ) -> ([f32; POINT_CNT * 3], [f32; POINT_CNT]) {
            let mut tmp_x = [0.0f32; POINT_CNT * 3];
            let mut tmp_y = [0.0f32; POINT_CNT];
            for i in 0..POINT_CNT {
                tmp_x[(i * 3 + 0) as usize] = 1.0f32;
                tmp_x[(i * 3 + 1) as usize] = x[i] as f32;
                tmp_x[(i * 3 + 2) as usize] = (x[i] as f32).powf(2.0f32);
                tmp_y[i] = (y[i] as f32).ln();
            }
            (tmp_x, tmp_y)
        }

        /// Solve ax = b.
        /// Any singular value smaller than `eps` is assumed to be zero.
        pub fn get_coefficents_from_8_matrix_data(
            &mut self,
            a: &[f32; POINT_CNT * 3],
            b: &[f32; POINT_CNT],
            eps: f32,
        ) {
            if a.len() == 3 * b.len() {
                type MatrixXx1f32 = SMatrix<f32, POINT_CNT, 1>;
                type MatrixXx3f32 = SMatrix<f32, POINT_CNT, 3>;
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
        }
    }
}
