use full_palette::{BLACK, BLUE, GREEN, RED};
use plotters::prelude::*;
use std::collections::HashMap;

use super::{atom::Atom, molecule::Molecule, recipient::Recipient}; // Plotting library

// Function to plot points

pub fn plot_recipient(recipient: &Recipient, output_path: &str, plot_grid: bool) {
    let root_area = BitMapBackend::new(
        output_path,
        (recipient.shape.0 as u32, recipient.shape.1 as u32),
    )
    .into_drawing_area();
    root_area.fill(&BLACK).unwrap();

    let mut chart = ChartBuilder::on(&root_area)
        .build_cartesian_2d(0.0..recipient.shape.0, -recipient.shape.1..0.0)
        .unwrap();

    chart.configure_mesh().draw().unwrap();

    for atom in recipient.contents.iter() {
        chart
            .draw_series(vec![Circle::new(
                (atom.position.0, -atom.position.1),
                atom.radius,
                ShapeStyle {
                    color: (atom.color).into(),
                    filled: true,
                    stroke_width: 1,
                },
            )])
            .unwrap();
    }
    if plot_grid {
        let (cols, rows) = recipient
            .grid
            .as_ref()
            .expect("Cant plot grid if no grid precent")
            .cells;
        let (width, height) = recipient.shape;
        for i in 0..=cols {
            let x = i as f32 * (width / cols as f32);
            chart
                .draw_series(LineSeries::new(
                    vec![(x, 0.0), (x, -height)],
                    WHITE.mix(0.8).stroke_width(1),
                ))
                .unwrap();
        }
        for j in 0..=rows {
            let y = -j as f32 * (height / rows as f32);
            chart
                .draw_series(LineSeries::new(
                    vec![(0.0, y), (width, y)],
                    WHITE.mix(0.8).stroke_width(1),
                ))
                .unwrap();
        }
    }
}

fn plot_paths(paths: Vec<(bool, Vec<(f64, f64)>)>, output_path: &str) {
    let root_area = BitMapBackend::new(output_path, (800, 600)).into_drawing_area();
    root_area.fill(&BLACK).unwrap();

    let mut chart = ChartBuilder::on(&root_area)
        .caption("Generated Path Points", ("Arial", 20))
        .build_cartesian_2d(0.0..1000.0, -1000.0..0.0) // Adjust the range as needed
        .unwrap();

    chart.configure_mesh().draw().unwrap();

    // Plot all the points
    for (i, path) in paths.iter().enumerate() {
        let color = {
            if (i % 2) == 0 {
                &RED
            } else {
                &BLUE
            }
        };

        let (closed, points) = path;
        chart
            .draw_series(PointSeries::of_element(
                points.iter().map(|p| (p.0, -p.1)),
                2,
                color,
                &|c, s, _| {
                    return EmptyElement::at(c)  // Position the points
                    + Circle::new((0,0), s, ShapeStyle {
                        color: color.to_rgba(),
                        filled: true,
                        stroke_width: 0,
                    });
                },
            ))
            .unwrap();
    }
}
