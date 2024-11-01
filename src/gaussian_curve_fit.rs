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

        // solve ax = b
        pub fn get_coefficents_from_8_matrix_data(
            &mut self,
            a: &[f32; POINT_CNT * 3],
            b: &[f32; POINT_CNT],
        ) {
            if a.len() == 3 * b.len() {
                type Matrix6x1f32 = SMatrix<f32, POINT_CNT, 1>;
                type Matrix6x3f32 = SMatrix<f32, POINT_CNT, 3>;
                let mb = Matrix6x1f32::from_row_slice(b);
                let ma = Matrix6x3f32::from_row_slice(a);
                // debug_rprintln!("{}", ma);

                let decomp = ma.svd(true, true);

                let x = decomp.solve(&mb, 1e-4);
                if let Ok(r) = x {
                    // debug_rprintln!("***{}", r);
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
