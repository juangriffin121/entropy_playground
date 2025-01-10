use std::fs;
use std::path::Path;
use svg_path_parser::parse_with_resolution;

pub fn extension(filepath: &str) -> Option<FileExtension> {
    match Path::new(filepath).extension() {
        Some(ext) if ext == "svg" => Some(FileExtension::SVG),
        _ => None,
    }
}

pub enum FileExtension {
    SVG,
}

pub fn load(filepath: &str) -> Option<Vec<(bool, Vec<(f64, f64)>)>> {
    let ext = extension(filepath);

    match ext {
        Some(e) => match e {
            FileExtension::SVG => Some(load_svg(filepath)),
        },
        _ => None,
    }
}

fn load_svg(filepath: &str) -> Vec<(bool, Vec<(f64, f64)>)> {
    // Read the SVG file
    let svg_data = fs::read_to_string(filepath).expect("Failed to read SVG file");

    // Extract the path data (the `d` attribute of the <path> tag) using a regex

    let re = regex::Regex::new(r#"<path[^>]*\sd="([^"]+)""#).expect("Failed to compile regex");
    let mut all_paths = Vec::new();

    for caps in re.captures_iter(&svg_data) {
        if let Some(d) = caps.get(1) {
            // Parse the path string

            let parsed_paths =
                parse_with_resolution(&d.as_str(), 1).collect::<Vec<(bool, Vec<(f64, f64)>)>>();
            all_paths.extend(parsed_paths);
        }
    }

    all_paths
}
