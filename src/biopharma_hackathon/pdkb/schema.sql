-- Derived views over the Parkinson's repurposing knowledge base.
--
-- The ingested tables are loaded by build.py; everything here is derived, so the
-- whole analytic layer is readable in one place and re-runnable after an agent
-- retunes score_weights.
--
--   drug_target_evidence  one row per (drug, target) with direction + evidence
--   drug_candidates       one row per drug: score components, score, verdict
--   undrugged_targets     PD targets with no approved drug, and screen coverage
--   pathway_evidence      per-pathway convergence of the two evidence routes

-- One row per (drug, target). The direction column is the gate: it says whether
-- this drug pushes the target the way that opposes the toxin insult, the way
-- that mimics it, or in a direction nobody can defend.
CREATE OR REPLACE VIEW drug_target_evidence AS
WITH pathway_rollup AS (
    SELECT protein_name AS gene,
           min(pathway_qvalue)                                     AS best_pathway_qvalue,
           count(DISTINCT pathway_index)                           AS n_pathways,
           count(DISTINCT CASE WHEN pathway_qvalue < 0.05
                               THEN pathway_index END)             AS n_enriched_pathways,
           bool_or(protein_is_pd_seed)                             AS is_pd_seed
    FROM pd_tree
    GROUP BY 1
),
toxin_rollup AS (
    SELECT gene_symbol                        AS gene,
           count(DISTINCT toxin)              AS n_toxins,
           sum(specificity_weight)            AS toxin_support,
           string_agg(DISTINCT toxin, '; ')   AS toxins
    FROM toxin_target
    GROUP BY 1
)
SELECT
    td.drug,
    td.gene,
    td.uniprot,
    td.protein,
    td.action_type,
    td.approval_year,
    -- protective / risk / ambiguous / unknown; see mechanism.py for the curation
    CASE
        WHEN md.gene IS NULL OR td.action_type IS NULL           THEN 'unknown'
        WHEN list_contains(str_split(md.protective_actions, '; '), td.action_type)
                                                                 THEN 'protective'
        WHEN list_contains(str_split(md.risk_actions, '; '), td.action_type)
                                                                 THEN 'risk'
        WHEN md.confidence = 'ambiguous'                         THEN 'ambiguous'
        ELSE 'unknown'
    END                                       AS direction,
    md.confidence                             AS direction_confidence,
    md.rationale                              AS direction_rationale,
    tr.n_toxins,
    tr.toxin_support,
    tr.toxins,
    pr.best_pathway_qvalue,
    pr.n_pathways,
    pr.n_enriched_pathways,
    coalesce(pr.is_pd_seed, FALSE)            AS target_is_pd_seed,
    sc.n_screen_hits,
    sc.best_screen_score
FROM target_drug td
LEFT JOIN mechanism_direction md ON md.gene = td.gene
LEFT JOIN toxin_rollup       tr ON tr.gene = td.gene
LEFT JOIN pathway_rollup     pr ON pr.gene = td.gene
LEFT JOIN (
    SELECT gene, count(*) AS n_screen_hits, max(score) AS best_screen_score
    FROM screen_hits GROUP BY 1
) sc ON sc.gene = td.gene
WHERE td.has_approved_drug;

-- Weights as a single row, so drug_candidates can cross join them. An agent can
-- UPDATE score_weights and re-select the view to see a different ranking.
CREATE OR REPLACE VIEW score_weight_row AS
SELECT
    max(CASE WHEN component = 'pathway'   THEN weight END) AS w_pathway,
    max(CASE WHEN component = 'toxin'     THEN weight END) AS w_toxin,
    max(CASE WHEN component = 'direction' THEN weight END) AS w_direction,
    max(CASE WHEN component = 'clinical'  THEN weight END) AS w_clinical,
    max(CASE WHEN component = 'biomarker' THEN weight END) AS w_biomarker,
    max(CASE WHEN component = 'cns'       THEN weight END) AS w_cns
FROM score_weights;

-- One row per drug. Every component is kept alongside the score so a ranking can
-- always be taken apart; the score is a convenience, not the evidence.
CREATE OR REPLACE VIEW drug_candidates AS
WITH per_drug AS (
    SELECT
        drug,
        count(DISTINCT gene)                                          AS n_targets,
        string_agg(DISTINCT gene, '; ')                               AS targets,
        min(approval_year)                                            AS approval_year,
        count(DISTINCT CASE WHEN direction = 'protective' THEN gene END) AS n_protective,
        count(DISTINCT CASE WHEN direction = 'risk'       THEN gene END) AS n_risk,
        count(DISTINCT CASE WHEN direction = 'ambiguous'  THEN gene END) AS n_ambiguous,
        min(best_pathway_qvalue)                                      AS best_pathway_qvalue,
        sum(n_enriched_pathways)                                      AS n_enriched_pathways,
        max(n_toxins)                                                 AS max_toxins_per_target,
        sum(toxin_support)                                            AS toxin_support,
        string_agg(DISTINCT toxins, '; ')                             AS toxins,
        bool_or(target_is_pd_seed)                                    AS hits_pd_seed_protein
    FROM drug_target_evidence
    GROUP BY 1
),
trial_rollup AS (
    SELECT drug,
           count(DISTINCT nct_id)                                            AS n_trials,
           count(DISTINCT CASE WHEN biomarker_measures THEN nct_id END)      AS n_biomarker_trials,
           count(DISTINCT CASE WHEN is_parkinsons THEN nct_id END)           AS n_pd_trials,
           count(DISTINCT CASE WHEN is_neuro      THEN nct_id END)           AS n_neuro_trials,
           count(DISTINCT CASE WHEN phase IN ('Phase 3', 'Phase 4')
                               THEN nct_id END)                              AS n_late_phase_trials
    FROM trials
    GROUP BY 1
),
raw AS (
    SELECT p.*,
           -- WHO INN stems: -mab is a monoclonal antibody, -cept an Fc fusion,
           -- -ase a therapeutic enzyme. All are large molecules that do not
           -- meaningfully cross the blood-brain barrier.
           CASE WHEN lower(p.drug) LIKE '%mab'
                  OR lower(p.drug) LIKE '%cept'
                  OR lower(p.drug) LIKE '%ase'  THEN 'biologic'
                ELSE 'small_molecule' END      AS modality,
           coalesce(t.n_trials, 0)            AS n_trials,
           coalesce(t.n_biomarker_trials, 0)  AS n_biomarker_trials,
           coalesce(t.n_pd_trials, 0)         AS n_pd_trials,
           coalesce(t.n_neuro_trials, 0)      AS n_neuro_trials,
           coalesce(t.n_late_phase_trials, 0) AS n_late_phase_trials
    FROM per_drug p
    LEFT JOIN trial_rollup t ON t.drug = p.drug
),
-- Components are normalised across the candidate set, so they are comparable
-- to each other but only meaningful relative to this cohort of 71 drugs.
scored AS (
    SELECT r.*,
           coalesce(
               -log10(nullif(r.best_pathway_qvalue, 0))
               / nullif(max(-log10(nullif(r.best_pathway_qvalue, 0))) OVER (), 0),
               0)                                                   AS pathway_component,
           coalesce(r.toxin_support / nullif(max(r.toxin_support) OVER (), 0), 0)
                                                                    AS toxin_component,
           -- maps [-1, +1] onto [0, 1]; 0.5 means no directional information
           ((r.n_protective - r.n_risk)::DOUBLE / nullif(r.n_targets, 0) + 1) / 2
                                                                    AS direction_component,
           ln(1 + r.n_trials) / nullif(max(ln(1 + r.n_trials)) OVER (), 0)
                                                                    AS clinical_component,
           coalesce(r.n_biomarker_trials::DOUBLE / nullif(r.n_trials, 0), 0)
                                                                    AS biomarker_component,
           -- Can the drug plausibly reach the brain? Without this, trial volume
           -- alone floats systemic blockbusters to the top: the anti-TNF
           -- antibodies have 100-469 trials each and zero in any neurological
           -- indication, yet PD is a CNS disease and IgG reaches the CSF at
           -- well under 1% of serum. Antibodies are not zeroed -- a peripheral
           -- immune mechanism is a live hypothesis -- only heavily discounted.
           CASE WHEN r.n_pd_trials > 0    THEN 1.00
                WHEN r.n_neuro_trials > 0 THEN 0.70
                WHEN r.modality = 'biologic' THEN 0.15
                ELSE 0.50 END                                       AS cns_component
    FROM raw r
)
SELECT
    s.drug,
    s.targets,
    s.n_targets,
    s.modality,
    s.toxins,
    s.approval_year,
    -- The verdict is deliberately not derived from the score. A drug that pushes
    -- a target the wrong way is not a weak candidate, it is the wrong hypothesis.
    CASE
        WHEN s.n_risk > 0 AND s.n_protective = 0 THEN 'mechanistic_risk'
        WHEN s.n_pd_trials >= 5                  THEN 'established_pd_therapy'
        WHEN s.n_protective > 0                  THEN 'candidate'
        WHEN s.n_ambiguous > 0                   THEN 'direction_ambiguous'
        ELSE 'direction_unknown'
    END                                          AS verdict,
    s.n_protective, s.n_risk, s.n_ambiguous,
    s.best_pathway_qvalue, s.n_enriched_pathways, s.hits_pd_seed_protein,
    s.toxin_support, s.max_toxins_per_target,
    s.n_trials, s.n_late_phase_trials, s.n_biomarker_trials,
    s.n_pd_trials, s.n_neuro_trials,
    round(s.pathway_component, 4)::DOUBLE   AS pathway_component,
    round(s.toxin_component, 4)::DOUBLE     AS toxin_component,
    round(s.direction_component, 4)::DOUBLE AS direction_component,
    round(s.clinical_component, 4)::DOUBLE  AS clinical_component,
    round(s.biomarker_component, 4)::DOUBLE AS biomarker_component,
    round(s.cns_component, 4)::DOUBLE       AS cns_component,
    round(
        w.w_pathway   * s.pathway_component
      + w.w_toxin     * s.toxin_component
      + w.w_direction * s.direction_component
      + w.w_clinical  * s.clinical_component
      + w.w_biomarker * s.biomarker_component
      + w.w_cns       * s.cns_component
    , 4)::DOUBLE                    AS priority_score
FROM scored s
CROSS JOIN score_weight_row w;

-- PD-implicated targets that no approved drug engages by direct mechanism. These
-- are where repurposing has nothing to offer and virtual screening might.
CREATE OR REPLACE VIEW undrugged_targets AS
SELECT
    td.gene,
    td.uniprot,
    td.protein,
    tt.n_toxins,
    tt.toxin_support,
    tt.toxins,
    pr.best_pathway_qvalue,
    pr.n_enriched_pathways,
    md.confidence  AS direction_confidence,
    md.protective_actions,
    coalesce(sc.n_screen_hits, 0) AS n_screen_hits,
    sc.best_screen_score
FROM (SELECT DISTINCT gene, uniprot, protein FROM target_drug WHERE NOT has_approved_drug) td
LEFT JOIN (
    SELECT gene_symbol AS gene, count(DISTINCT toxin) AS n_toxins,
           sum(specificity_weight) AS toxin_support, string_agg(DISTINCT toxin, '; ') AS toxins
    FROM toxin_target GROUP BY 1
) tt ON tt.gene = td.gene
LEFT JOIN (
    SELECT protein_name AS gene, min(pathway_qvalue) AS best_pathway_qvalue,
           count(DISTINCT CASE WHEN pathway_qvalue < 0.05 THEN pathway_index END)
               AS n_enriched_pathways
    FROM pd_tree GROUP BY 1
) pr ON pr.gene = td.gene
LEFT JOIN mechanism_direction md ON md.gene = td.gene
LEFT JOIN (
    SELECT gene, count(*) AS n_screen_hits, max(score) AS best_screen_score
    FROM screen_hits GROUP BY 1
) sc ON sc.gene = td.gene;

-- Where the disease-association route and the toxicological route agree.
CREATE OR REPLACE VIEW pathway_evidence AS
SELECT
    e.pathway_index,
    e.pathway_name,
    e.pathway_size,
    e.qvalue           AS pathway_qvalue,
    e.fold_enrichment,
    e.seed_overlap,
    coalesce(tp.n_toxins, 0)       AS n_toxins,
    coalesce(tp.n_toxin_targets, 0) AS n_toxin_targets,
    coalesce(tp.toxin_support, 0)  AS toxin_support,
    tp.toxins,
    tp.toxin_targets,
    (e.qvalue < 0.05 AND tp.pathway_index IS NOT NULL) AS is_convergent
FROM pathway_enrichment e
LEFT JOIN toxin_pathway tp USING (pathway_index);
