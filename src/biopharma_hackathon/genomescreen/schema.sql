-- GenomeScreen relational schema.
--
-- Grain:
--   targets           one row per UniProt accession
--   pockets           one row per screen-result directory (target x fragment x pocket)
--   pocket_structures one row per refined receptor conformation + its docking grid
--   molecules         one row per distinct (catalog id, SMILES) pair
--   hits              one row per (pocket, molecule) DrugCLIP hit

CREATE TABLE targets (
    uniprot_acc      VARCHAR PRIMARY KEY,
    n_pockets        INTEGER NOT NULL,
    n_hits           BIGINT  NOT NULL,
    best_score       DOUBLE
);

CREATE TABLE pockets (
    pocket_key       VARCHAR PRIMARY KEY,  -- the screen_results/ directory name
    uniprot_acc      VARCHAR NOT NULL,
    af_entry         INTEGER NOT NULL,     -- the F<n> in the AlphaFold DB id
    af_model_version INTEGER NOT NULL,
    fragment_idx     INTEGER NOT NULL,     -- AF2 model fragment, for proteins split by length
    pocket_kind      VARCHAR NOT NULL,     -- 'template' (from an aligned holo PDB) | 'detected'
    pocket_idx       INTEGER NOT NULL,
    -- derived by the rollup pass at the end of the build
    n_structures     INTEGER NOT NULL DEFAULT 0,
    n_hits           INTEGER NOT NULL DEFAULT 0,
    best_score       DOUBLE,
    dir_path         VARCHAR NOT NULL
);

CREATE TABLE pocket_structures (
    pocket_key       VARCHAR NOT NULL,
    structure_idx    INTEGER NOT NULL,
    template_pdb_id  VARCHAR,             -- set only for pocket_kind = 'template'
    structure_path   VARCHAR NOT NULL,    -- gzipped PDB of the refined complex
    grid_path        VARCHAR,             -- NULL for the 3 structures shipped without a grid
    grid_center_x    DOUBLE,
    grid_center_y    DOUBLE,
    grid_center_z    DOUBLE,
    grid_file        VARCHAR,
    receptor_file    VARCHAR,
    PRIMARY KEY (pocket_key, structure_idx)
);

CREATE TABLE molecules (
    mol_id           INTEGER PRIMARY KEY,
    oid              VARCHAR NOT NULL,    -- vendor id exactly as shipped
    smiles           VARCHAR NOT NULL,
    source           VARCHAR NOT NULL,    -- 'ZINC' | 'Enamine REAL' | 'Enamine PV' | 'unknown'
    catalog_id       VARCHAR NOT NULL,    -- oid without the protomer/tautomer suffix
    protomer_idx     INTEGER,             -- Enamine only
    tautomer_idx     INTEGER,             -- Enamine only
    UNIQUE (oid, smiles)
);

CREATE TABLE hits (
    pocket_key       VARCHAR NOT NULL,
    mol_id           INTEGER NOT NULL,
    score            DOUBLE  NOT NULL,    -- DrugCLIP similarity
    rank_in_pocket   INTEGER NOT NULL,    -- 1 = best-scoring hit for this pocket
    source_index     BIGINT  NOT NULL     -- leader.csv 'Name': row index in the screened library
);

CREATE INDEX idx_pockets_acc      ON pockets (uniprot_acc);
CREATE INDEX idx_structures_pdb   ON pocket_structures (template_pdb_id);
CREATE INDEX idx_molecules_oid    ON molecules (oid);
CREATE INDEX idx_molecules_cat    ON molecules (catalog_id);
CREATE INDEX idx_hits_pocket      ON hits (pocket_key);
CREATE INDEX idx_hits_mol         ON hits (mol_id);

-- Flat view for ad-hoc querying: every hit with its target and molecule attached.
CREATE VIEW hit_details AS
SELECT
    p.uniprot_acc,
    h.pocket_key,
    p.pocket_kind,
    p.fragment_idx,
    p.pocket_idx,
    m.oid,
    m.catalog_id,
    m.source,
    m.smiles,
    h.score,
    h.rank_in_pocket
FROM hits AS h
JOIN pockets   AS p USING (pocket_key)
JOIN molecules AS m USING (mol_id);

-- Best score per (target, molecule), collapsing pockets and fragments.
CREATE VIEW target_molecule_best AS
SELECT
    p.uniprot_acc,
    m.catalog_id,
    m.source,
    max(h.score)          AS best_score,
    count(*)              AS n_hits,
    count(DISTINCT h.pocket_key) AS n_pockets,
    arg_max(m.smiles, h.score)   AS best_smiles
FROM hits AS h
JOIN pockets   AS p USING (pocket_key)
JOIN molecules AS m USING (mol_id)
GROUP BY 1, 2, 3;
