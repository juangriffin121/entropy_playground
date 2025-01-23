use std::collections::HashSet;

use super::{
    atom::{DistributionJSON, ElementsJSON},
    bond::BondDefinitions,
    molecule::MoleculesJSON,
    recipient::BlueprintAtomIndex,
};

pub struct JSONChecker {}

impl JSONChecker {
    pub fn all_checks(
        elements: &ElementsJSON,
        distribution: &Option<DistributionJSON>,
        molecules: &Option<MoleculesJSON>,
        bonds: &Option<BondDefinitions>,
    ) {
        if let Some(distribution) = distribution {
            match JSONChecker::check_distribution(&elements, &distribution) {
                Ok(_) => println!("Distribution is consistent with ElementsJSON!"),
                Err(err) => eprintln!("Consistency check failed: {}", err),
            }
        }

        if let Some(molecules) = molecules {
            match JSONChecker::check_molecules_bonds(molecules) {
                Ok(_) => println!("All molecule bonds are valid!"),
                Err(err) => eprintln!("Molecule bond check failed: {}", err),
            }

            match JSONChecker::check_molecules_atoms(elements, molecules) {
                Ok(_) => println!("All molecule atoms are consistent with ElementsJSON!"),
                Err(err) => eprintln!("Molecule atom check failed: {}", err),
            }
        }

        if let Some(bonds) = bonds {
            match JSONChecker::check_bond_definitions(&elements, &bonds) {
                Ok(_) => println!("Bond definitions are consistent with ElementsJSON!"),
                Err(err) => eprintln!("Bond definitions consistency check failed: {}", err),
            }
        }
    }
    pub fn check_distribution(
        elements: &ElementsJSON,
        distribution: &DistributionJSON,
    ) -> Result<(), String> {
        for atom in distribution.map.keys() {
            if !elements.data.contains_key(atom) {
                return Err(format!(
                    "Atom '{}' in distribution is not found in elements.json",
                    atom
                ));
            }
        }
        Ok(())
    }

    pub fn check_molecules_bonds(molecules: &MoleculesJSON) -> Result<(), String> {
        for (molecule_name, molecule) in &molecules.blueprints {
            let mut seen_bonds = HashSet::new();
            let mut bonded_atoms: HashSet<BlueprintAtomIndex> = HashSet::new();
            for (idx, &(atom1_idx, atom2_idx)) in molecule.bonds.iter().enumerate() {
                if atom1_idx == atom2_idx {
                    return Err(format!(
                        "Invalid bond in molecule {}: bond from atom {} to itself",
                        molecule_name, atom1_idx.0
                    ));
                }

                if atom1_idx.0 >= molecule.atoms.len() || atom2_idx.0 >= molecule.atoms.len() {
                    return Err(format!(
                        "Invalid bond in molecule {}: 
                            bond references atom index out of range: ({}, {})",
                        molecule_name, atom1_idx.0, atom2_idx.0
                    ));
                }
                let ordered_bond = (
                    BlueprintAtomIndex(atom1_idx.0.min(atom2_idx.0)),
                    BlueprintAtomIndex(atom1_idx.0.max(atom2_idx.0)),
                );
                if !seen_bonds.insert(ordered_bond) {
                    return Err(format!(
                        "Duplicate bond found in molecule {}: {:?}",
                        molecule_name, ordered_bond
                    ));
                }
                bonded_atoms.insert(atom1_idx);
                bonded_atoms.insert(atom2_idx);
            }
            for (idx, _) in molecule.atoms.iter().enumerate() {
                if !bonded_atoms.contains(&BlueprintAtomIndex(idx)) {
                    return Err(format!(
                        "Invalid molecule {}: atom {} is not part of any bond",
                        molecule_name, idx
                    ));
                }
            }
        }
        Ok(())
    }

    pub fn check_molecules_atoms(
        elements: &ElementsJSON,
        molecules: &MoleculesJSON,
    ) -> Result<(), String> {
        for (name, blueprint) in &molecules.blueprints {
            for atom in &blueprint.atoms {
                if !elements.data.contains_key(atom) {
                    return Err(format!(
                        "Atom '{}' in molecule '{}' is not found in elements.json",
                        atom, name
                    ));
                }
            }
        }
        Ok(())
    }
    pub fn check_bond_definitions(
        elements: &ElementsJSON,
        bonds: &BondDefinitions,
    ) -> Result<(), String> {
        let tuple_map = bonds.to_tuple_map();
        for ((element1, element2), _) in tuple_map.iter() {
            if !elements.data.contains_key(element1) || !elements.data.contains_key(element2) {
                return Err(format!(
                    "Bond definition contains unknown element(s): {}-{}",
                    element1, element2
                ));
            }
        }
        Ok(())
    }
}
