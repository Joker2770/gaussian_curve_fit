# gaussian_curve_fit

[![Rust](https://github.com/Joker2770/gaussian_curve_fit/actions/workflows/rust.yml/badge.svg)](https://github.com/Joker2770/gaussian_curve_fit/actions/workflows/rust.yml)

A `no_std` and no `alloc` library for gaussian curve coefficents calculation.

## example

```rust
    let mut gaussian_coes = GaussianCoefficents2D::default();
    let xdata = [
        -8.0f32, -6.0f32, -4.0f32, -2.0f32, 0.0f32, 2.0f32, 4.0f32, 6.0f32,
    ];
    let ydata = [
        6.7f32, 10.6f32, 13.5f32, 15.7f32, 16.6f32, 15.4f32, 14.2f32, 10.3f32,
    ];
    let (x_arr, y_arr) = GaussianCoefficents2D::get_matrix_data_from_8_points(&xdata, &ydata);
    let _ = gaussian_coes.get_coefficents_from_8_matrix_data(&x_arr, &y_arr, 1e-4);

    assert!((gaussian_coes.value(-8.0f32).unwrap_or_default() - 6.7f32).abs() < 1.0f32);
    assert!((gaussian_coes.value(-4.0f32).unwrap_or_default() - 13.5f32).abs() < 1.0f32);
    assert!((gaussian_coes.value(0.0f32).unwrap_or_default() - 16.6f32).abs() < 1.0f32);
    assert!((gaussian_coes.value(4.0f32).unwrap_or_default() - 14.2f32).abs() < 1.0f32);

```
