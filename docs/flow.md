flowchart TD
  subgraph Classification
    A1[DefaultClassification]  
    A2[UserDictClassification.from_json]
    A1 --> B[.classify(resname,atom)]
    A2 --> B
  end

  subgraph Topology
    B --> C[read_top(resname, .rtf/.str)]
    C --> D[radii[resname] = {atom: [radius, code]}]
  end

  subgraph Trajectory
    E[run_defect.py] -->|load radii| B
    E --> F[mda.Universe(top,traj)]
    F --> G[select_atoms('resname…')]
    G --> H[PackingDefect2Sequential]
  end

  subgraph Per-Frame Loop
    H --> I{frame ts}
    I --> J1[apply_pbc]
    I --> J2[calculate_bounding_box]
    I --> J3[initialize_grid]
    I --> J4[select leaflets]
    J4 --> J5[populate_grid_with_atoms]  
    J5 --> I
  end

  subgraph Conclude & Output
    I --> K[aggregate grids & zlims]
    K --> L[initialize_empty_defect_universe]
    L --> M[place defect atoms at (x,y,z)]
    M --> N[write_combined_gro(protein_atoms,defect_atoms)]
  end

