`%>%` <- magrittr::`%>%`
library(patchwork)

# Utility functions ------------------------------------------------------------

read_m8 <- function(file) {
  cols <- readr::cols(
    query = readr::col_character(),
    subject = readr::col_character(),
    pid = readr::col_double(),
    ali_length = readr::col_double(),
    mismatches = readr::col_double(),
    gap_openings = readr::col_double(),
    query_start = readr::col_double(),
    query_end = readr::col_double(),
    subject_start = readr::col_double(),
    subject_end = readr::col_double(),
    evalue = readr::col_double(),
    bitscore = readr::col_double()
  )
  out <- readr::read_tsv(
    file,
    col_names = names(cols$cols),
    col_types = cols,
    progress = FALSE
  )
  out
}

read_tadb_meta <- function(path) {
  cols <- readr::cols(
    type = readr::col_character(),
    validation = readr::col_character(),
    id = readr::col_character(),
    related_mge = readr::col_character(),
    organism = readr::col_character(),
    seqid = readr::col_character(),
    toxin = readr::col_character(),
    antitoxin = readr::col_character(),
    family = readr::col_character()
  )
  raw <- readr::read_tsv(
    path,
    col_names = names(cols$cols),
    col_types = cols$cols,
    skip = 1,
    progress = FALSE
  )
  out <- raw %>%
    tidyr::pivot_longer(
      cols = c(toxin, antitoxin),
      names_to = "component",
      values_to = "gene"
    ) %>%
    tidyr::separate_wider_delim(
      cols = family,
      delim = "/",
      names = c("family", "domain")
    ) %>%
    dplyr::mutate(product = stringr::str_extract(gene, "(?<=\\().*?(?=\\))"))
}

read_gff <- function(file) {
  cols <- readr::cols(
    seqid = readr::col_character(),
    source = readr::col_character(),
    feature = readr::col_character(),
    start = readr::col_double(),
    end = readr::col_double(),
    score = readr::col_double(),
    strand = readr::col_character(),
    frame = readr::col_character(),
    attribute = readr::col_character()
  )
  out <- readr::read_tsv(
    file,
    col_names = names(cols$cols),
    col_types = cols,
    progress = FALSE,
    comment = "#"
  )
  out
}

# Filter plasmids --------------------------------------------------------------

fs::dir_create("output")

gkab859_sup_path <- "data/gkab859_supplementary_data.xlsx"

gkab859_sup_tbl_01 <- gkab859_sup_path %>%
  readxl::read_xlsx(sheet = 2, progress = FALSE)

plasmid_representatives <- gkab859_sup_tbl_01 %>%
  janitor::clean_names() %>%
  dplyr::filter(derep == TRUE) %>%
  dplyr::select(acc) %>%
  dplyr::pull()

plasmid_representatives_nover <- plasmid_representatives %>%
  stringr::str_remove_all("\\..*")

# CRISPR-Cas on plasmids -------------------------------------------------------

# Cas operons

gkab859_sup_path <- "data/gkab859_supplementary_data.xlsx"

gkab859_sup_tbl_05 <- gkab859_sup_path %>%
  readxl::read_xlsx(sheet = 6, progress = FALSE)

cas_on_plasmids <- gkab859_sup_tbl_05 %>%
  janitor::clean_names() %>%
  dplyr::select(acc, prediction)

# All plasmids with Cas operons

plasmids_w_cas <- cas_on_plasmids %>%
  dplyr::distinct(acc) %>%
  dplyr::pull(acc)

plasmids_w_cas_nover <- plasmids_w_cas %>%
  stringr::str_remove_all("\\..*")

# All plasmids with Cas operons, excluding Type IV

plasmids_w_cas_no_iv <- cas_on_plasmids %>%
  dplyr::filter(!stringr::str_detect(prediction, "^IV-")) %>%
  dplyr::distinct(acc) %>%
  dplyr::pull(acc)

plasmids_w_cas_no_iv_nover <- plasmids_w_cas_no_iv %>%
  stringr::str_remove_all("\\..*")

# All plasmids with Cas operons, only Type IV

plasmids_w_cas_only_iv <- cas_on_plasmids %>%
  dplyr::filter(stringr::str_detect(prediction, "^IV-")) %>%
  dplyr::distinct(acc) %>%
  dplyr::pull(acc)

plasmids_w_cas_only_iv_nover <- plasmids_w_cas_only_iv %>%
  stringr::str_remove_all("\\..*")

cas_on_plasmids_summary <- cas_on_plasmids %>%
  dplyr::summarise(n = dplyr::n(), .by = prediction) %>%
  dplyr::arrange(dplyr::desc(n))

cas_on_plasmids_summary %>%
  dplyr::mutate(
    type = dplyr::if_else(stringr::str_detect(prediction, "IV-"), "IV", "Other")
  ) %>%
  dplyr::summarise(n = sum(n), .by = type) %>%
  dplyr::mutate(p = n / sum(n))

cas_on_plasmids_summary %>%
  dplyr::rename(`Cas operon type` = prediction) %>%
  dplyr::arrange(`Cas operon type`) %>%
  writexl::write_xlsx("output/cas_counts.xlsx")

# CRISPR operons

gkab859_sup_tbl_03 <- gkab859_sup_path %>%
  readxl::read_xlsx(sheet = 4, progress = FALSE)

crispr_on_plasmids <- gkab859_sup_tbl_03 %>%
  janitor::clean_names() %>%
  dplyr::select(acc, subtype)

crispr_on_plasmids_summary <- crispr_on_plasmids %>%
  dplyr::summarise(n = dplyr::n(), .by = subtype) %>%
  dplyr::arrange(dplyr::desc(n))

crispr_on_plasmids_summary %>%
  dplyr::mutate(
    subtype = dplyr::if_else(subtype == "NA", "Unknown", subtype)
  ) %>%
  dplyr::arrange(subtype) %>%
  dplyr::rename(`CRISPR array type` = subtype) %>%
  writexl::write_xlsx("output/crispr_counts.xlsx")

# Plots

plot_01 <- cas_on_plasmids_summary %>%
  dplyr::mutate(
    Category = dplyr::case_when(
      stringr::str_detect(prediction, "IV") ~ "Type IV",
      .default = "Other"
    )
  ) %>%
  ggplot2::ggplot(ggplot2::aes(x = prediction, y = n, fill = Category)) +
  ggplot2::geom_col(colour = "black", linewidth = 0.24) +
  ggplot2::scale_fill_manual(values = c("#e5e5e9", "#e9c54e")) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.1))) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 45,
      hjust = 1,
      size = 5,
      colour = "black"
    ),
    legend.position = "top",
    text = ggplot2::element_text(size = 9, colour = "black"),
    axis.text = ggplot2::element_text(size = 9, colour = "black"),
    line = ggplot2::element_line(linewidth = 72.27 / 96 * 0.5),
    rect = ggplot2::element_rect(linewidth = 72.27 / 96 * 0.5),
    panel.border = ggplot2::element_rect(linewidth = 72.27 / 96 * 0.5),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.minor.x = ggplot2::element_blank(),
    panel.grid.minor.y = ggplot2::element_blank()
  ) +
  ggplot2::labs(
    y = "Count",
    x = "Cas operon type",
    fill = ""
  ) +
  ggplot2::coord_cartesian(clip = "off")

plot_01

plot_02_x_order <- crispr_on_plasmids_summary %>%
  dplyr::filter(subtype != "NA") %>%
  dplyr::arrange(subtype) %>%
  dplyr::pull(subtype) %>%
  c(., "Unknown")

plot_02 <- crispr_on_plasmids_summary %>%
  dplyr::mutate(
    subtype = dplyr::if_else(subtype == "NA", "Unknown", subtype)
  ) %>%
  dplyr::mutate(
    Category = dplyr::case_when(
      stringr::str_detect(subtype, "IV") ~ "Type IV",
      .default = "Other"
    )
  ) %>%
  ggplot2::ggplot(ggplot2::aes(
    x = factor(subtype, levels = plot_02_x_order),
    y = n,
    fill = Category
  )) +
  ggplot2::geom_col(colour = "black", linewidth = 0.24) +
  ggplot2::scale_fill_manual(values = c("#e5e5e9", "#e9c54e")) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.1))) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 5),
    legend.position = "top",
    text = ggplot2::element_text(size = 9, colour = "black"),
    axis.text = ggplot2::element_text(size = 9, colour = "black"),
    line = ggplot2::element_line(linewidth = 72.27 / 96 * 0.5),
    rect = ggplot2::element_rect(linewidth = 72.27 / 96 * 0.5),
    panel.border = ggplot2::element_rect(linewidth = 72.27 / 96 * 0.5),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.minor.x = ggplot2::element_blank(),
    panel.grid.minor.y = ggplot2::element_blank()
  ) +
  ggplot2::labs(
    y = "Count",
    x = "CRISPR array type"
  ) +
  ggplot2::coord_cartesian(clip = "off")

plot_02

fig <- plot_01 + plot_02 + patchwork::plot_layout()

fig

fig %>%
  ggplot2::ggsave(
    filename = "output/crispr_cas_counts.pdf",
    width = 182.4,
    height = 65,
    units = "mm"
  )


# Targeting network ------------------------------------------------------------

targeting_network_path <- "data/gkab859_supplementary_spacers.tsv"

targeting_network <- targeting_network_path %>%
  readr::read_tsv()

plasmid_targeting_network <- targeting_network %>%
  janitor::clean_names() %>%
  dplyr::mutate(
    source = stringr::str_remove_all(spacer, "_[0-9]*@.*")
  ) %>%
  dplyr::filter(
    source_type == "Plasmid" &
      target_type == "Plasmid"
  ) %>%
  dplyr::select(
    source,
    target
  )

plasmid_targeting_network_filt <- plasmid_targeting_network %>%
  dplyr::filter(
    source %in%
      plasmid_representatives_nover &
      target %in% plasmid_representatives_nover &
      source != target
  )

plasmid_targeting_network_filt_no_iv <- plasmid_targeting_network %>%
  dplyr::filter(
    source %in% plasmids_w_cas_no_iv_nover & source != target
  )

plasmid_targeting_network_filt_only_iv <- plasmid_targeting_network %>%
  dplyr::filter(
    source %in% plasmids_w_cas_only_iv_nover & source != target
  )

targeted_plasmids <- plasmid_targeting_network_filt %>%
  dplyr::distinct(target) %>%
  dplyr::pull()

targeted_plasmids_no_iv <- plasmid_targeting_network_filt_no_iv %>%
  dplyr::distinct(target) %>%
  dplyr::pull()

targeted_plasmids_only_iv <- plasmid_targeting_network_filt_only_iv %>%
  dplyr::distinct(target) %>%
  dplyr::pull()

# Get plasmid ids --------------------------------------------------------------

plasmids_fna_path <- "data/plasmid_seq.fna"

plasmids_fna <- Biostrings::readDNAStringSet(plasmids_fna_path)

plasmid_ids <- plasmids_fna %>%
  names() %>%
  stringr::str_extract("^\\S+") %>%
  dplyr::as_tibble() %>%
  dplyr::rename(plasmid_id = value) %>%
  dplyr::mutate(
    plasmid_id_nover = stringr::str_remove_all(plasmid_id, "\\..*")
  )

# Filter TA hits ---------------------------------------------------------------

tadb_meta_path <- "data/tadb_meta.tsv"
tadb_meta_raw <- read_tadb_meta(tadb_meta_path)

tadb_meta <- tadb_meta_raw %>%
  dplyr::mutate(
    query = dplyr::if_else(
      component == "toxin",
      id %>% stringr::str_remove("TA[0]*") %>% paste0("T", .),
      id %>% stringr::str_remove("TA[0]*") %>% paste0("AT", .)
    )
  ) %>%
  dplyr::select(query, id, type, family, domain, product) %>%
  dplyr::rename(ta_type = type)

aa_m8_path <- "data/ta_blast/plasmid_faa.m8"
na_m8_path <- "data/ta_blast/plasmid_fna.m8"

aa_m8_raw <- read_m8(aa_m8_path)
na_m8_raw <- read_m8(na_m8_path)

gff_path <- "data/plasmid_seq.gff"

plasmid_genes <- gff_path %>%
  read_gff() %>%
  dplyr::arrange(seqid, start, end) %>%
  dplyr::mutate(relative_position = dplyr::row_number(), .by = seqid) %>%
  dplyr::mutate(subject = paste0(seqid, "_", relative_position))

query_na_path <- "data/tadb.fna"
query_na <- Biostrings::readDNAStringSet(query_na_path)

query_names <- query_na %>% names() %>% stringr::str_extract("^\\S+")
query_lengths <- query_na %>% Biostrings::width()
query_lengths <- tibble::tibble(
  query = query_names,
  query_length = query_lengths
)

na_m8 <- na_m8_raw %>%
  dplyr::select(-c(mismatches, gap_openings)) %>%
  dplyr::left_join(query_lengths, by = dplyr::join_by(query)) %>%
  dplyr::mutate(query_ali_length = query_end + 1 - query_start) %>%
  dplyr::mutate(query_coverage = query_ali_length / query_length) %>%
  dplyr::filter(query_coverage >= 0.7) %>%
  # Filter overlapping hits for best scoring
  dplyr::arrange(subject, dplyr::desc(bitscore)) %>%
  dplyr::filter(
    !duplicated(cumsum(c(0, diff(subject_start) > 0))),
    .by = subject
  ) %>%
  dplyr::mutate(start = subject_start, end = subject_end)

query_aa_path <- "data/tadb.faa"
query_aa <- Biostrings::readAAStringSet(query_aa_path)

query_names <- query_aa %>% names() %>% stringr::str_extract("^\\S+")
query_lengths <- query_aa %>% Biostrings::width()
query_lengths <- tibble::tibble(
  query = query_names,
  query_length = query_lengths
)

subject_aa_path <- "data/plasmid_seq.faa"
subject_aa <- Biostrings::readAAStringSet(subject_aa_path)

subject_names <- subject_aa %>% names() %>% stringr::str_extract("^\\S+")
subject_lengths <- subject_aa %>% Biostrings::width()
subject_lengths <- tibble::tibble(
  subject = subject_names,
  subject_length = subject_lengths
)

aa_m8 <- aa_m8_raw %>%
  dplyr::left_join(query_lengths, by = dplyr::join_by(query)) %>%
  dplyr::left_join(subject_lengths, by = dplyr::join_by(subject)) %>%
  dplyr::mutate(query_ali_length = query_end + 1 - query_start) %>%
  dplyr::mutate(query_coverage = query_ali_length / query_length) %>%
  dplyr::mutate(subject_ali_length = subject_end + 1 - subject_start) %>%
  dplyr::mutate(subject_coverage = subject_ali_length / subject_length) %>%
  dplyr::filter(query_coverage >= 0.7 & subject_coverage >= 0.7) %>%
  dplyr::left_join(
    plasmid_genes %>% dplyr::select(subject, start, end),
    by = dplyr::join_by(subject)
  )

hits <- dplyr::bind_rows(aa_m8, na_m8) %>%
  dplyr::mutate(seqid = stringr::str_remove(subject, "_([^_]*)$")) %>%
  dplyr::select(-subject)

plasmid_genes_new <- na_m8 %>%
  dplyr::select(subject, subject_start, subject_end) %>%
  dplyr::mutate(type = "ncrna") %>%
  dplyr::rename(seqid = subject, start = subject_start, end = subject_end) %>%
  dplyr::arrange(seqid, start, end) %>%
  dplyr::mutate(
    subject = paste0(seqid, "_ncrna_", dplyr::row_number()),
    .by = seqid
  ) %>%
  dplyr::bind_rows(plasmid_genes) %>%
  dplyr::arrange(seqid, start, end) %>%
  dplyr::mutate(relative_position = dplyr::row_number(), .by = seqid)

hits_new <- plasmid_genes_new %>%
  dplyr::left_join(
    hits,
    by = dplyr::join_by(seqid, start, end),
    relationship = "many-to-many"
  )

hits_best <- hits_new %>%
  dplyr::filter(!is.na(query)) %>%
  dplyr::left_join(tadb_meta, by = dplyr::join_by(query)) %>%
  dplyr::select(
    seqid,
    subject,
    start,
    end,
    relative_position,
    query,
    bitscore,
    query_coverage,
    subject_coverage,
    ta_type,
    family
  ) %>%
  dplyr::slice_max(
    order_by = bitscore,
    n = 1,
    by = c(subject, family),
    with_ties = FALSE
  ) %>%
  dplyr::arrange(seqid, family, relative_position) %>%
  dplyr::mutate(
    family_cluster = cumsum(c(1, abs(diff(relative_position)) > 0 + 1)),
    .by = c(seqid)
  ) %>%
  dplyr::mutate(
    family_score = sum(
      bitscore * mean(c(query_coverage, subject_coverage), na.rm = TRUE),
      na.rm = TRUE
    ),
    .by = c(seqid, family, family_cluster)
  ) %>%
  dplyr::slice_max(
    order_by = family_score,
    n = 1,
    by = c(subject),
    with_ties = FALSE
  )

ta_candidates <- hits_best %>%
  dplyr::mutate(
    component = dplyr::if_else(
      stringr::str_detect(query, "^T"),
      "toxin",
      "antitoxin"
    )
  ) %>%
  dplyr::mutate(
    candidate = dplyr::if_else(
      "toxin" %in% component & "antitoxin" %in% component,
      TRUE,
      FALSE
    ),
    .by = c(seqid, family_cluster, family)
  ) %>%
  dplyr::filter(candidate == TRUE)

plasmids_w_ta <- ta_candidates %>%
  dplyr::distinct(seqid) %>%
  dplyr::pull()

plasmids_w_ta_nover <- plasmids_w_ta %>%
  stringr::str_remove_all("\\..*")

# Plot TA hits -----------------------------------------------------------------

ta_candidates %>%
  dplyr::distinct(seqid, family) %>%
  dplyr::summarise(n = dplyr::n(), .by = family) %>%
  dplyr::arrange(dplyr::desc(n))

ta_count <- ta_candidates %>%
  dplyr::filter(seqid %in% plasmid_representatives) %>%
  dplyr::arrange(seqid, family, relative_position) %>%
  dplyr::mutate(
    ta_system = sort(unique(family)) %>% paste0(collapse = "/"),
    ta_type = sort(unique(ta_type)) %>% paste0(collapse = "/"),
    .by = c(seqid, family_cluster, family)
  ) %>%
  dplyr::distinct(seqid, ta_type, ta_system) %>%
  dplyr::summarise(n = dplyr::n(), .by = c(ta_type, ta_system)) %>%
  dplyr::arrange(dplyr::desc(n)) %>%
  dplyr::filter(ta_type != "")

ta_count %>%
  dplyr::filter(ta_system == "parDE")

order <- ta_count %>% dplyr::pull(ta_system)

psks <- c("hok-sok", "ccdAB", "pemIK", "mazEF", "abiEi-abiEii", "toxIN")

plot_03 <- ta_count %>%
  dplyr::mutate(
    Category = dplyr::case_when(
      ta_system == "parDE" ~ "ParDE",
      ta_system %in% psks ~ "Known PSK system",
      .default = "Other"
    )
  ) %>%
  ggplot2::ggplot(ggplot2::aes(x = factor(ta_system), y = n, fill = Category)) +
  ggplot2::geom_col(colour = "black", linewidth = 0.24) +
  ggplot2::scale_fill_manual(values = c("#b35807", "#e5e5e9", "#7fc97f")) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.05))) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    legend.position = "top"
  ) +
  ggplot2::labs(
    x = "Toxin-antitoxin family",
    y = "Number"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 6,
      colour = "black"
    ),
    # axis.title.y = ggplot2::element_blank(),
    # axis.text.y = ggplot2::element_blank(),
    # axis.ticks.y = ggplot2::element_blank(),
    text = ggplot2::element_text(size = 9, colour = "black"),
    axis.text = ggplot2::element_text(size = 9, colour = "black"),
    line = ggplot2::element_line(linewidth = 72.27 / 96 * 0.5),
    rect = ggplot2::element_rect(linewidth = 72.27 / 96 * 0.5),
    panel.border = ggplot2::element_rect(linewidth = 72.27 / 96 * 0.5),
    # panel.grid.major.y = ggplot2::element_blank(),
    panel.grid.minor.x = ggplot2::element_blank(),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.minor.y = ggplot2::element_blank(),
    strip.background = ggplot2::element_blank(),
    legend.position = "top"
  ) +
  ggplot2::facet_grid(. ~ ta_type, scales = "free_x", space = "free")

plot_03

ta_proportion_targeted <- ta_candidates %>%
  dplyr::arrange(seqid, family, relative_position) %>%
  dplyr::mutate(
    ta_system = sort(unique(family)) %>% paste0(collapse = "/"),
    ta_type = sort(unique(ta_type)) %>% paste0(collapse = "/"),
    .by = c(seqid, family_cluster, family)
  ) %>%
  dplyr::distinct(seqid, ta_type, ta_system) %>%
  dplyr::mutate(seqid_nover = stringr::str_remove_all(seqid, "\\..*")) %>%
  dplyr::mutate(
    is_targeted = dplyr::if_else(
      seqid_nover %in% targeted_plasmids,
      TRUE,
      FALSE
    )
  ) %>%
  dplyr::summarise(n = dplyr::n(), .by = c(ta_type, ta_system, is_targeted)) %>%
  dplyr::mutate(p = n / sum(n), .by = c(ta_type, ta_system)) %>%
  dplyr::filter(ta_type != "") %>%
  tidyr::complete(ta_system, is_targeted, fill = list(n = 0, p = 0)) %>%
  dplyr::mutate(ta_type = max(ta_type, na.rm = TRUE), .by = ta_system) %>%
  dplyr::filter(is_targeted == TRUE)

ta_proportion_targeted

plot_04 <- ta_proportion_targeted %>%
  dplyr::mutate(
    Category = dplyr::case_when(
      ta_system == "parDE" ~ "ParDE",
      ta_system %in% psks ~ "Known PSK system",
      .default = "Other"
    )
  ) %>%
  ggplot2::ggplot(ggplot2::aes(x = factor(ta_system), y = p, fill = Category)) +
  ggplot2::geom_col(colour = "black", linewidth = 0.24) +
  ggplot2::scale_fill_manual(values = c("#b35807", "#e5e5e9", "#7fc97f")) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(mult = c(0, 0)),
    limits = c(0, 1),
    labels = scales::label_percent()
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 45,
      hjust = 1,
      size = 9,
      colour = "black"
    ),
    # axis.title.y = ggplot2::element_blank(),
    # axis.text.y = ggplot2::element_blank(),
    # axis.ticks.y = ggplot2::element_blank(),
    text = ggplot2::element_text(size = 9, colour = "black"),
    axis.text = ggplot2::element_text(size = 9, colour = "black"),
    line = ggplot2::element_line(linewidth = 72.27 / 96 * 0.5),
    rect = ggplot2::element_rect(linewidth = 72.27 / 96 * 0.5),
    panel.border = ggplot2::element_rect(linewidth = 72.27 / 96 * 0.5),
    # panel.grid.major.y = ggplot2::element_blank(),
    panel.grid.minor.x = ggplot2::element_blank(),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.minor.y = ggplot2::element_blank(),
    strip.background = ggplot2::element_blank(),
    legend.position = "top"
  ) +
  ggplot2::labs(
    y = "Proportion",
    x = "Toxin-antitoxin family"
  ) +
  ggplot2::facet_grid(. ~ ta_type, scales = "free_x", space = "free")

plot_04

ta_proportion_targeted_no_type_iv <- ta_candidates %>%
  dplyr::filter(seqid %in% plasmid_representatives) %>%
  dplyr::arrange(seqid, family, relative_position) %>%
  dplyr::mutate(
    ta_system = sort(unique(family)) %>% paste0(collapse = "/"),
    ta_type = sort(unique(ta_type)) %>% paste0(collapse = "/"),
    .by = c(seqid, family_cluster, family)
  ) %>%
  dplyr::distinct(seqid, ta_type, ta_system) %>%
  dplyr::mutate(seqid_nover = stringr::str_remove_all(seqid, "\\..*")) %>%
  dplyr::mutate(
    is_targeted = dplyr::if_else(
      seqid_nover %in% targeted_plasmids_no_iv,
      TRUE,
      FALSE
    )
  ) %>%
  dplyr::summarise(n = dplyr::n(), .by = c(ta_type, ta_system, is_targeted)) %>%
  dplyr::mutate(p = n / sum(n), .by = c(ta_type, ta_system)) %>%
  dplyr::filter(ta_type != "") %>%
  tidyr::complete(ta_system, is_targeted, fill = list(n = 0, p = 0)) %>%
  dplyr::mutate(ta_type = max(ta_type, na.rm = TRUE), .by = ta_system) %>%
  dplyr::filter(is_targeted == TRUE)

plot_05 <- ta_proportion_targeted_no_type_iv %>%
  dplyr::mutate(
    Category = dplyr::case_when(
      ta_system == "parDE" ~ "ParDE",
      ta_system %in% psks ~ "Known PSK system",
      .default = "Other"
    )
  ) %>%
  ggplot2::ggplot(ggplot2::aes(x = factor(ta_system), y = p, fill = Category)) +
  ggplot2::geom_col(colour = "black", linewidth = 0.24) +
  ggplot2::scale_fill_manual(values = c("#b35807", "#e5e5e9", "#7fc97f")) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(mult = c(0, 0)),
    limits = c(0, 1),
    labels = scales::label_percent()
  ) +
  ggplot2::labs(
    x = "Toxin-antitoxin family",
    y = "Proportion"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 6,
      colour = "black"
    ),
    # axis.title.y = ggplot2::element_blank(),
    # axis.text.y = ggplot2::element_blank(),
    # axis.ticks.y = ggplot2::element_blank(),
    text = ggplot2::element_text(size = 9, colour = "black"),
    axis.text = ggplot2::element_text(size = 9, colour = "black"),
    line = ggplot2::element_line(linewidth = 72.27 / 96 * 0.5),
    rect = ggplot2::element_rect(linewidth = 72.27 / 96 * 0.5),
    panel.border = ggplot2::element_rect(linewidth = 72.27 / 96 * 0.5),
    # panel.grid.major.y = ggplot2::element_blank(),
    panel.grid.minor.x = ggplot2::element_blank(),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.minor.y = ggplot2::element_blank(),
    strip.background = ggplot2::element_blank(),
    legend.position = "top"
  ) +
  ggplot2::facet_grid(. ~ ta_type, scales = "free_x", space = "free")

plot_05

ta_proportion_targeted_only_type_iv <- ta_candidates %>%
  dplyr::filter(seqid %in% plasmid_representatives) %>%
  dplyr::arrange(seqid, family, relative_position) %>%
  dplyr::mutate(
    ta_system = sort(unique(family)) %>% paste0(collapse = "/"),
    ta_type = sort(unique(ta_type)) %>% paste0(collapse = "/"),
    .by = c(seqid, family_cluster, family)
  ) %>%
  dplyr::distinct(seqid, ta_type, ta_system) %>%
  dplyr::mutate(seqid_nover = stringr::str_remove_all(seqid, "\\..*")) %>%
  dplyr::mutate(
    is_targeted = dplyr::if_else(
      seqid_nover %in% targeted_plasmids_only_iv,
      TRUE,
      FALSE
    )
  ) %>%
  dplyr::summarise(n = dplyr::n(), .by = c(ta_type, ta_system, is_targeted)) %>%
  dplyr::mutate(p = n / sum(n), .by = c(ta_type, ta_system)) %>%
  dplyr::filter(ta_type != "") %>%
  tidyr::complete(ta_system, is_targeted, fill = list(n = 0, p = 0)) %>%
  dplyr::mutate(ta_type = max(ta_type, na.rm = TRUE), .by = ta_system) %>%
  dplyr::filter(is_targeted == TRUE)

plot_06 <- ta_proportion_targeted_only_type_iv %>%
  dplyr::mutate(
    Category = dplyr::case_when(
      ta_system == "parDE" ~ "ParDE",
      ta_system %in% psks ~ "Known PSK system",
      .default = "Other"
    )
  ) %>%
  ggplot2::ggplot(ggplot2::aes(x = factor(ta_system), y = p, fill = Category)) +
  ggplot2::geom_col(colour = "black", linewidth = 0.24) +
  ggplot2::scale_fill_manual(values = c("#b35807", "#e5e5e9", "#7fc97f")) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(mult = c(0, 0)),
    limits = c(0, 1),
    labels = scales::label_percent()
  ) +
  ggplot2::labs(
    x = "Toxin-antitoxin family",
    y = "Proportion"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 6,
      colour = "black"
    ),
    # axis.title.y = ggplot2::element_blank(),
    # axis.text.y = ggplot2::element_blank(),
    # axis.ticks.y = ggplot2::element_blank(),
    text = ggplot2::element_text(size = 9, colour = "black"),
    axis.text = ggplot2::element_text(size = 9, colour = "black"),
    line = ggplot2::element_line(linewidth = 72.27 / 96 * 0.5),
    rect = ggplot2::element_rect(linewidth = 72.27 / 96 * 0.5),
    panel.border = ggplot2::element_rect(linewidth = 72.27 / 96 * 0.5),
    # panel.grid.major.y = ggplot2::element_blank(),
    panel.grid.minor.x = ggplot2::element_blank(),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.minor.y = ggplot2::element_blank(),
    strip.background = ggplot2::element_blank(),
    legend.position = "top"
  ) +
  ggplot2::facet_grid(. ~ ta_type, scales = "free_x", space = "free")

plot_06

layout <- "
A
B
C
"

fig <- plot_03 +
  plot_06 +
  plot_05 +
  patchwork::plot_layout(design = layout)

fig

fig %>%
  ggplot2::ggsave(
    filename = "output/ta_counts.pdf",
    width = 182.4,
    height = 250,
    units = "mm"
  )

ta_count %>%
  dplyr::left_join(
    ta_proportion_targeted %>%
      dplyr::select(-c(is_targeted, n)) %>%
      dplyr::rename(`Proportion targeted (including Type IV)` = p)
  ) %>%
  dplyr::left_join(
    ta_proportion_targeted_no_type_iv %>%
      dplyr::select(-c(is_targeted, n)) %>%
      dplyr::rename(`Proportion targeted (excluding Type IV)` = p)
  ) %>%
  dplyr::rename(`TA type` = ta_type, `TA family` = ta_system) %>%
  dplyr::arrange(`TA type`, `TA family`) %>%
  writexl::write_xlsx("output/ta_counts.xlsx")

# Statistical analysis ---------------------------------------------------------

# All

plasmid_features_all <- plasmid_ids %>%
  dplyr::mutate(
    is_targeted = dplyr::if_else(
      plasmid_id_nover %in% targeted_plasmids,
      TRUE,
      FALSE
    ),
    has_ta = dplyr::if_else(
      plasmid_id %in% plasmids_w_ta,
      TRUE,
      FALSE
    )
  )

contingency_table_all <- table(
  plasmid_features_all$has_ta,
  plasmid_features_all$is_targeted
)

rownames(contingency_table_all) <- c(FALSE, TRUE)
colnames(contingency_table_all) <- c(FALSE, TRUE)

chi_test_all <- chisq.test(contingency_table_all)
chi_test_all

residuals_matrix_all <- chi_test_all$residuals
residuals_matrix_all

pvals_matrix_all <- 2 * pnorm(-abs(residuals_matrix_all))
pvals_matrix_all

expected_all <- chi_test_all$expected

expected_se_all <- sqrt(expected_all)

N_all <- sum(contingency_table_all)

phi_all <- sqrt(chi_test_all$statistic / N_all)
phi_all

observed_counts_all <- contingency_table_all %>%
  as.data.frame() %>%
  dplyr::as_tibble() %>%
  dplyr::rename(has_ta = Var1, targeted = Var2, count = Freq) %>%
  dplyr::mutate(type = "Observed")

expected_counts_all <- observed_counts_all %>%
  dplyr::mutate(count = as.vector(expected_all)) %>%
  dplyr::mutate(se = as.vector(expected_se_all)) %>%
  dplyr::mutate(type = "Expected")

counts_all <- dplyr::bind_rows(observed_counts_all, expected_counts_all)

o_all <- observed_counts_all %>%
  dplyr::rename(observed_count = count) %>%
  dplyr::select(-type)

e_all <- expected_counts_all %>%
  dplyr::rename(expected_count = count) %>%
  dplyr::select(-type)

e_vs_o_all <- o_all %>%
  dplyr::left_join(e_all, by = dplyr::join_by(has_ta, targeted)) %>%
  dplyr::mutate(
    label = paste0(
      scales::comma(observed_count),
      "\n(",
      scales::comma(round(expected_count)),
      " ±",
      round(se),
      ")"
    )
  ) %>%
  dplyr::mutate(
    direction = dplyr::case_when(
      dplyr::between(observed_count, expected_count - se, expected_count + se) ~
        "No difference",
      observed_count > expected_count ~ "More than expected",
      observed_count < expected_count ~ "Less than expected",
    )
  ) %>%
  dplyr::mutate(group = "All")

plot_01 <- e_vs_o_all %>%
  ggplot2::ggplot(ggplot2::aes(x = targeted, y = has_ta)) +
  ggplot2::geom_tile(ggplot2::aes(fill = direction), colour = "black") +
  ggplot2::geom_text(
    ggplot2::aes(label = label, colour = direction),
    size = 2.5
  ) +
  ggplot2::scale_fill_manual(values = c("#f6ceca", "#c5e5fb", "white")) +
  ggplot2::scale_colour_manual(values = c("#9c251c", "#01488d", "black")) +
  ggplot2::labs(x = "Plasmid is targeted", y = "Plasmid has TA") +
  ggplot2::scale_x_discrete(expand = c(0, 0)) +
  ggplot2::scale_y_discrete(expand = c(0, 0)) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.text.y = ggplot2::element_text(angle = 90, hjust = 0.5),
    text = ggplot2::element_text(size = 9, colour = "black"),
    axis.text = ggplot2::element_text(size = 9, colour = "black"),
    line = ggplot2::element_line(linewidth = 72.27 / 96 * 0.5),
    rect = ggplot2::element_rect(linewidth = 72.27 / 96 * 0.5),
    panel.border = ggplot2::element_rect(linewidth = 72.27 / 96 * 0.5),
    axis.ticks = ggplot2::element_blank(),
    legend.position = "NA"
  ) +
  ggplot2::coord_cartesian(clip = "off")

plot_01

plot_02 <- observed_counts_all %>%
  dplyr::summarise(count = sum(count), .by = has_ta) %>%
  ggplot2::ggplot(ggplot2::aes(x = count, y = has_ta)) +
  ggplot2::geom_col(
    fill = "#e5e5e9",
    colour = "black",
    linewidth = 72.27 / 96 * 0.5
  ) +
  ggplot2::scale_x_continuous(
    expand = ggplot2::expansion(mult = c(0, 0)),
    limits = c(0, 19000),
    n.breaks = 4
  ) +
  ggplot2::labs(x = "Count", y = NULL) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 45,
      hjust = 1,
      size = 9,
      colour = "black"
    ),
    axis.title.y = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    text = ggplot2::element_text(size = 9, colour = "black"),
    axis.text = ggplot2::element_text(size = 9, colour = "black"),
    line = ggplot2::element_line(linewidth = 72.27 / 96 * 0.5),
    rect = ggplot2::element_rect(linewidth = 72.27 / 96 * 0.5),
    panel.border = ggplot2::element_rect(linewidth = 72.27 / 96 * 0.5),
    panel.grid.major.y = ggplot2::element_blank(),
    panel.grid.minor.x = ggplot2::element_blank(),
    panel.grid.minor.y = ggplot2::element_blank()
  ) +
  ggplot2::coord_cartesian(clip = "off")

plot_02

plot_03 <- observed_counts_all %>%
  dplyr::summarise(count = sum(count), .by = targeted) %>%
  ggplot2::ggplot(ggplot2::aes(y = count, x = targeted)) +
  ggplot2::geom_col(
    fill = "#e5e5e9",
    colour = "black",
    linewidth = 72.27 / 96 * 0.5
  ) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(mult = c(0, 0)),
    limits = c(0, 19000),
    n.breaks = 4
  ) +
  ggplot2::labs(y = "Count", y = NULL) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.title.x = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank(),
    text = ggplot2::element_text(size = 9, colour = "black"),
    axis.text = ggplot2::element_text(size = 9, colour = "black"),
    line = ggplot2::element_line(linewidth = 72.27 / 96 * 0.5),
    rect = ggplot2::element_rect(linewidth = 72.27 / 96 * 0.5),
    panel.border = ggplot2::element_rect(linewidth = 72.27 / 96 * 0.5),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.minor.x = ggplot2::element_blank(),
    panel.grid.minor.y = ggplot2::element_blank()
  ) +
  ggplot2::coord_cartesian(clip = "off")

plot_03

layout <- "
CC#
AAB
AAB
"

plot_04 <- patchwork::wrap_plots(plot_01, plot_02, plot_03, design = layout)

plot_04

plot_04 %>%
  ggplot2::ggsave(
    filename = "output/ct_all.pdf",
    width = 89,
    height = 89,
    units = "mm"
  )

# Excluding Type IV

plasmid_features_no_iv <- plasmid_ids %>%
  dplyr::mutate(
    is_targeted = dplyr::if_else(
      plasmid_id_nover %in% targeted_plasmids_no_iv,
      TRUE,
      FALSE
    ),
    has_ta = dplyr::if_else(
      plasmid_id %in% plasmids_w_ta,
      TRUE,
      FALSE
    )
  )

contingency_table_no_iv <- table(
  plasmid_features_no_iv$has_ta,
  plasmid_features_no_iv$is_targeted
)

rownames(contingency_table_no_iv) <- c(FALSE, TRUE)
colnames(contingency_table_no_iv) <- c(FALSE, TRUE)

chi_test_no_iv <- chisq.test(contingency_table_no_iv)
chi_test_no_iv

residuals_matrix_no_iv <- chi_test_no_iv$residuals
residuals_matrix_no_iv

pvals_matrix_no_iv <- 2 * pnorm(-abs(residuals_matrix_no_iv))
pvals_matrix_no_iv

expected_no_iv <- chi_test_no_iv$expected

expected_se_no_iv <- sqrt(expected_no_iv)

N_no_iv <- sum(contingency_table_no_iv)

phi_no_iv <- sqrt(chi_test_no_iv$statistic / N_no_iv)
phi_no_iv

observed_counts_no_iv <- contingency_table_no_iv %>%
  as.data.frame() %>%
  dplyr::as_tibble() %>%
  dplyr::rename(has_ta = Var1, targeted = Var2, count = Freq) %>%
  dplyr::mutate(type = "Observed")

expected_counts_no_iv <- observed_counts_no_iv %>%
  dplyr::mutate(count = as.vector(expected_no_iv)) %>%
  dplyr::mutate(se = as.vector(expected_se_no_iv)) %>%
  dplyr::mutate(type = "Expected")

counts_no_iv <- dplyr::bind_rows(observed_counts_no_iv, expected_counts_no_iv)

o_no_iv <- observed_counts_no_iv %>%
  dplyr::rename(observed_count = count) %>%
  dplyr::select(-type)

e_no_iv <- expected_counts_no_iv %>%
  dplyr::rename(expected_count = count) %>%
  dplyr::select(-type)

e_vs_o_no_iv <- o_no_iv %>%
  dplyr::left_join(e_no_iv, by = dplyr::join_by(has_ta, targeted)) %>%
  dplyr::mutate(
    label = paste0(
      scales::comma(observed_count),
      "\n(",
      scales::comma(round(expected_count)),
      " ±",
      round(se),
      ")"
    )
  ) %>%
  dplyr::mutate(
    direction = dplyr::case_when(
      dplyr::between(observed_count, expected_count - se, expected_count + se) ~
        "No difference",
      observed_count > expected_count ~ "More than expected",
      observed_count < expected_count ~ "Less than expected",
    )
  ) %>%
  dplyr::mutate(group = "Non-type IV")

plot_01 <- e_vs_o_no_iv %>%
  ggplot2::ggplot(ggplot2::aes(x = targeted, y = has_ta)) +
  ggplot2::geom_tile(ggplot2::aes(fill = direction), colour = "black") +
  ggplot2::geom_text(
    ggplot2::aes(label = label, colour = direction),
    size = 2.5
  ) +
  ggplot2::scale_fill_manual(values = c("#f6ceca", "#c5e5fb", "white")) +
  ggplot2::scale_colour_manual(values = c("#9c251c", "#01488d", "black")) +
  ggplot2::labs(x = "Plasmid is targeted", y = "Plasmid has TA") +
  ggplot2::scale_x_discrete(expand = c(0, 0)) +
  ggplot2::scale_y_discrete(expand = c(0, 0)) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.text.y = ggplot2::element_text(angle = 90, hjust = 0.5),
    text = ggplot2::element_text(size = 9, colour = "black"),
    axis.text = ggplot2::element_text(size = 9, colour = "black"),
    line = ggplot2::element_line(linewidth = 72.27 / 96 * 0.5),
    rect = ggplot2::element_rect(linewidth = 72.27 / 96 * 0.5),
    panel.border = ggplot2::element_rect(linewidth = 72.27 / 96 * 0.5),
    axis.ticks = ggplot2::element_blank(),
    legend.position = "NA"
  ) +
  ggplot2::coord_cartesian(clip = "off")

plot_01

plot_02 <- observed_counts_no_iv %>%
  dplyr::summarise(count = sum(count), .by = has_ta) %>%
  ggplot2::ggplot(ggplot2::aes(x = count, y = has_ta)) +
  ggplot2::geom_col(
    fill = "#e5e5e9",
    colour = "black",
    linewidth = 72.27 / 96 * 0.5
  ) +
  ggplot2::scale_x_continuous(
    expand = ggplot2::expansion(mult = c(0, 0)),
    limits = c(0, 19000),
    n.breaks = 4
  ) +
  ggplot2::labs(x = "Count", y = NULL) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 45,
      hjust = 1,
      size = 9,
      colour = "black"
    ),
    axis.title.y = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    text = ggplot2::element_text(size = 9, colour = "black"),
    axis.text = ggplot2::element_text(size = 9, colour = "black"),
    line = ggplot2::element_line(linewidth = 72.27 / 96 * 0.5),
    rect = ggplot2::element_rect(linewidth = 72.27 / 96 * 0.5),
    panel.border = ggplot2::element_rect(linewidth = 72.27 / 96 * 0.5),
    panel.grid.major.y = ggplot2::element_blank(),
    panel.grid.minor.x = ggplot2::element_blank(),
    panel.grid.minor.y = ggplot2::element_blank()
  ) +
  ggplot2::coord_cartesian(clip = "off")

plot_02

plot_03 <- observed_counts_no_iv %>%
  dplyr::summarise(count = sum(count), .by = targeted) %>%
  ggplot2::ggplot(ggplot2::aes(y = count, x = targeted)) +
  ggplot2::geom_col(
    fill = "#e5e5e9",
    colour = "black",
    linewidth = 72.27 / 96 * 0.5
  ) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(mult = c(0, 0)),
    limits = c(0, 19000),
    n.breaks = 4
  ) +
  ggplot2::labs(y = "Count", y = NULL) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.title.x = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank(),
    text = ggplot2::element_text(size = 9, colour = "black"),
    axis.text = ggplot2::element_text(size = 9, colour = "black"),
    line = ggplot2::element_line(linewidth = 72.27 / 96 * 0.5),
    rect = ggplot2::element_rect(linewidth = 72.27 / 96 * 0.5),
    panel.border = ggplot2::element_rect(linewidth = 72.27 / 96 * 0.5),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.minor.x = ggplot2::element_blank(),
    panel.grid.minor.y = ggplot2::element_blank()
  ) +
  ggplot2::coord_cartesian(clip = "off")

plot_03

layout <- "
CC#
AAB
AAB
"

plot_04 <- patchwork::wrap_plots(plot_01, plot_02, plot_03, design = layout)

plot_04

plot_04 %>%
  ggplot2::ggsave(
    filename = "output/ct_excluding_type_iv.pdf",
    width = 89,
    height = 89,
    units = "mm"
  )

# Only Type IV

plasmid_features_only_iv <- plasmid_ids %>%
  dplyr::mutate(
    is_targeted = dplyr::if_else(
      plasmid_id_nover %in% targeted_plasmids_only_iv,
      TRUE,
      FALSE
    ),
    has_ta = dplyr::if_else(
      plasmid_id %in% plasmids_w_ta,
      TRUE,
      FALSE
    )
  )

contingency_table_only_iv <- table(
  plasmid_features_only_iv$has_ta,
  plasmid_features_only_iv$is_targeted
)

rownames(contingency_table_only_iv) <- c(FALSE, TRUE)
colnames(contingency_table_only_iv) <- c(FALSE, TRUE)

chi_test_only_iv <- chisq.test(contingency_table_only_iv)
chi_test_only_iv

residuals_matrix_only_iv <- chi_test_only_iv$residuals
residuals_matrix_only_iv

pvals_matrix_only_iv <- 2 * pnorm(-abs(residuals_matrix_only_iv))
pvals_matrix_only_iv

expected_only_iv <- chi_test_only_iv$expected

expected_se_only_iv <- sqrt(expected_only_iv)

N_only_iv <- sum(contingency_table_only_iv)

phi_only_iv <- sqrt(chi_test_only_iv$statistic / N_only_iv)
phi_only_iv

observed_counts_only_iv <- contingency_table_only_iv %>%
  as.data.frame() %>%
  dplyr::as_tibble() %>%
  dplyr::rename(has_ta = Var1, targeted = Var2, count = Freq) %>%
  dplyr::mutate(type = "Observed")

expected_counts_only_iv <- observed_counts_only_iv %>%
  dplyr::mutate(count = as.vector(expected_only_iv)) %>%
  dplyr::mutate(se = as.vector(expected_se_only_iv)) %>%
  dplyr::mutate(type = "Expected")

counts_only_iv <- dplyr::bind_rows(
  observed_counts_only_iv,
  expected_counts_only_iv
)

o_only_iv <- observed_counts_only_iv %>%
  dplyr::rename(observed_count = count) %>%
  dplyr::select(-type)

e_only_iv <- expected_counts_only_iv %>%
  dplyr::rename(expected_count = count) %>%
  dplyr::select(-type)

e_vs_o_only_iv <- o_only_iv %>%
  dplyr::left_join(e_only_iv, by = dplyr::join_by(has_ta, targeted)) %>%
  dplyr::mutate(
    label = paste0(
      scales::comma(observed_count),
      "\n(",
      scales::comma(round(expected_count)),
      " ±",
      round(se),
      ")"
    )
  ) %>%
  dplyr::mutate(
    direction = dplyr::case_when(
      dplyr::between(observed_count, expected_count - se, expected_count + se) ~
        "No difference",
      observed_count > expected_count ~ "More than expected",
      observed_count < expected_count ~ "Less than expected",
    )
  ) %>%
  dplyr::mutate(group = "Type IV")

plot_01 <- e_vs_o_only_iv %>%
  ggplot2::ggplot(ggplot2::aes(x = targeted, y = has_ta)) +
  ggplot2::geom_tile(ggplot2::aes(fill = direction), colour = "black") +
  ggplot2::geom_text(
    ggplot2::aes(label = label, colour = direction),
    size = 2.5
  ) +
  ggplot2::scale_fill_manual(values = c("#f6ceca", "#c5e5fb", "white")) +
  ggplot2::scale_colour_manual(values = c("#9c251c", "#01488d", "black")) +
  ggplot2::labs(x = "Plasmid is targeted", y = "Plasmid has TA") +
  ggplot2::scale_x_discrete(expand = c(0, 0)) +
  ggplot2::scale_y_discrete(expand = c(0, 0)) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.text.y = ggplot2::element_text(angle = 90, hjust = 0.5),
    text = ggplot2::element_text(size = 9, colour = "black"),
    axis.text = ggplot2::element_text(size = 9, colour = "black"),
    line = ggplot2::element_line(linewidth = 72.27 / 96 * 0.5),
    rect = ggplot2::element_rect(linewidth = 72.27 / 96 * 0.5),
    panel.border = ggplot2::element_rect(linewidth = 72.27 / 96 * 0.5),
    axis.ticks = ggplot2::element_blank(),
    legend.position = "NA"
  ) +
  ggplot2::coord_cartesian(clip = "off")

plot_01

plot_02 <- observed_counts_only_iv %>%
  dplyr::summarise(count = sum(count), .by = has_ta) %>%
  ggplot2::ggplot(ggplot2::aes(x = count, y = has_ta)) +
  ggplot2::geom_col(
    fill = "#e5e5e9",
    colour = "black",
    linewidth = 72.27 / 96 * 0.5
  ) +
  ggplot2::scale_x_continuous(
    expand = ggplot2::expansion(mult = c(0, 0)),
    limits = c(0, 19000),
    n.breaks = 4
  ) +
  ggplot2::labs(x = "Count", y = NULL) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 45,
      hjust = 1,
      size = 9,
      colour = "black"
    ),
    axis.title.y = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    text = ggplot2::element_text(size = 9, colour = "black"),
    axis.text = ggplot2::element_text(size = 9, colour = "black"),
    line = ggplot2::element_line(linewidth = 72.27 / 96 * 0.5),
    rect = ggplot2::element_rect(linewidth = 72.27 / 96 * 0.5),
    panel.border = ggplot2::element_rect(linewidth = 72.27 / 96 * 0.5),
    panel.grid.major.y = ggplot2::element_blank(),
    panel.grid.minor.x = ggplot2::element_blank(),
    panel.grid.minor.y = ggplot2::element_blank()
  ) +
  ggplot2::coord_cartesian(clip = "off")

plot_02

plot_03 <- observed_counts_only_iv %>%
  dplyr::summarise(count = sum(count), .by = targeted) %>%
  ggplot2::ggplot(ggplot2::aes(y = count, x = targeted)) +
  ggplot2::geom_col(
    fill = "#e5e5e9",
    colour = "black",
    linewidth = 72.27 / 96 * 0.5
  ) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(mult = c(0, 0)),
    limits = c(0, 19000),
    n.breaks = 4
  ) +
  ggplot2::labs(y = "Count", y = NULL) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.title.x = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank(),
    text = ggplot2::element_text(size = 9, colour = "black"),
    axis.text = ggplot2::element_text(size = 9, colour = "black"),
    line = ggplot2::element_line(linewidth = 72.27 / 96 * 0.5),
    rect = ggplot2::element_rect(linewidth = 72.27 / 96 * 0.5),
    panel.border = ggplot2::element_rect(linewidth = 72.27 / 96 * 0.5),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.minor.x = ggplot2::element_blank(),
    panel.grid.minor.y = ggplot2::element_blank()
  ) +
  ggplot2::coord_cartesian(clip = "off")

plot_03

layout <- "
CC#
AAB
AAB
"

plot_04 <- patchwork::wrap_plots(plot_01, plot_02, plot_03, design = layout)

plot_04

plot_04 %>%
  ggplot2::ggsave(
    filename = "output/ct_only_type_iv.pdf",
    width = 89,
    height = 89,
    units = "mm"
  )

# Proportion of targeted / non-targeted plasmids encoding TA -------------------

stats <- dplyr::bind_rows(
  as.data.frame(as.table(residuals_matrix_all)) %>%
    dplyr::rename(has_ta = Var1, targeted = Var2, residual = Freq) %>%
    dplyr::mutate(group = "All"),
  as.data.frame(as.table(residuals_matrix_only_iv)) %>%
    dplyr::rename(has_ta = Var1, targeted = Var2, residual = Freq) %>%
    dplyr::mutate(group = "Type IV"),
  as.data.frame(as.table(residuals_matrix_no_iv)) %>%
    dplyr::rename(has_ta = Var1, targeted = Var2, residual = Freq) %>%
    dplyr::mutate(group = "Non-type IV")
) %>%
  dplyr::mutate(p_value = 2 * pnorm(-abs(residual))) %>%
  dplyr::mutate(
    significance = dplyr::case_when(
      p_value < 0.001 ~ "< 0.001",
      p_value < 0.01 ~ "< 0.01",
      p_value < 0.05 ~ "< 0.05",
      .default = "N.S."
    )
  )

stats

stats_summary <- dplyr::bind_rows(
  e_vs_o_all,
  e_vs_o_only_iv,
  e_vs_o_no_iv
) %>%
  dplyr::left_join(stats, by = dplyr::join_by(has_ta, targeted, group)) %>%
  dplyr::select(group, dplyr::everything()) %>%
  dplyr::rename(standard_error = se) %>%
  dplyr::select(-label)

stats_summary

stats_summary %>%
  writexl::write_xlsx("output/expected_vs_observed_stats.xlsx")

x_order <- c("All", "Type IV", "Non-type IV")

stats_summary_long <- stats_summary %>%
  tidyr::pivot_longer(
    cols = c(observed_count, expected_count),
    names_to = "category",
    values_to = "n"
  ) %>%
  dplyr::mutate(
    category = dplyr::if_else(
      category == "observed_count",
      "Observed",
      "Expected"
    )
  ) %>%
  dplyr::mutate(
    standard_error = dplyr::if_else(category == "Observed", NA, standard_error)
  ) %>%
  dplyr::mutate(
    direction = dplyr::if_else(category == "Observed", direction, "Expected")
  ) %>%
  dplyr::mutate(
    lower = n - standard_error,
    upper = n + standard_error
  ) %>%
  dplyr::mutate(
    quadrant = dplyr::case_when(
      has_ta == TRUE & targeted == TRUE ~ "TA+ & targeted",
      has_ta == FALSE & targeted == TRUE ~ "TA- & targeted",
      has_ta == TRUE & targeted == FALSE ~ "TA+ & not targeted",
      has_ta == FALSE & targeted == FALSE ~ "TA- & not targeted"
    )
  )

format_p <- function(p) {
  ifelse(
    p > 0.001,
    formatC(p, format = "f", digits = 4), # fixed decimal, 4 places
    gsub(
      "e\\-0*(\\d+)",
      " × 10<sup>−\\1</sup>",
      formatC(p, format = "e", digits = 2)
    )
  )
}

annotations <- stats_summary_long %>%
  dplyr::mutate(label_y = max(n) * 1.2, .by = quadrant) %>%
  dplyr::filter(category == "Observed") %>%
  dplyr::mutate(
    formatted_p = format_p(p_value),
    label = paste0(
      "Residual = ",
      round(residual, 2),
      "<br>",
      "<i>p</i> = ",
      formatted_p
    )
  ) %>%
  dplyr::select(group, quadrant, label, label_y)

dodge_width <- 0.9
bar_offset <- dodge_width / 2

lines <- stats_summary_long %>%
  dplyr::mutate(y = max(n) * 1.1, .by = quadrant) %>%
  dplyr::mutate(x_pos = match(group, x_order)) %>%
  dplyr::reframe(
    x_start = unique(x_pos) - bar_offset / 2,
    x_end = unique(x_pos) + bar_offset / 2,
    .by = c(group, quadrant, y)
  )

ticks <- lines %>%
  tidyr::pivot_longer(
    cols = c(x_start, x_end),
    names_to = "end",
    values_to = "x"
  ) %>%
  dplyr::mutate(
    yend = y - 0.02 * y
  )

plot_01 <- stats_summary_long %>%
  ggplot2::ggplot(ggplot2::aes(
    x = factor(group, x_order),
    y = n,
    fill = direction
  )) +
  ggplot2::geom_col(
    colour = "black",
    linewidth = 0.24,
    position = ggplot2::position_dodge(0.9)
  ) +
  ggplot2::scale_fill_manual(values = c("#e5e5e9", "#dc6465", "#5496ce")) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = lower, ymax = upper),
    position = ggplot2::position_dodge(0.9),
    width = 0.2,
    linewidth = 0.24
  ) +
  ggtext::geom_richtext(
    data = annotations,
    ggplot2::aes(x = group, y = label_y, label = label),
    inherit.aes = FALSE,
    fill = NA,
    label.color = NA,
    size = 5 / ggplot2::.pt,
    hjust = 0.5
  ) +
  ggplot2::geom_segment(
    data = lines,
    ggplot2::aes(x = x_start, xend = x_end, y = y, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.24,
    color = "black"
  ) +
  ggplot2::geom_segment(
    data = ticks,
    ggplot2::aes(x = x, xend = x, y = y, yend = yend),
    inherit.aes = FALSE,
    linewidth = 0.24,
    color = "black"
  ) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(mult = c(0, 0.1))
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 0,
      hjust = 0.5,
      size = 9,
      colour = "black"
    ),
    legend.position = "top",
    text = ggplot2::element_text(size = 9, colour = "black"),
    axis.text = ggplot2::element_text(size = 9, colour = "black"),
    line = ggplot2::element_line(linewidth = 72.27 / 96 * 0.5),
    rect = ggplot2::element_rect(linewidth = 72.27 / 96 * 0.5),
    panel.border = ggplot2::element_rect(linewidth = 72.27 / 96 * 0.5),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.minor.x = ggplot2::element_blank(),
    panel.grid.minor.y = ggplot2::element_blank()
  ) +
  ggplot2::labs(
    y = "Count",
    x = "Type of targeting CRISPR-Cas"
  ) +
  ggplot2::coord_cartesian(clip = "off") +
  ggplot2::facet_wrap(~quadrant, scales = "free")

plot_01

plot_01 %>%
  ggplot2::ggsave(
    filename = "output/expected_vs_observed.pdf",
    width = 182.4,
    height = 150,
    units = "mm"
  )

stats_summary_long_TA <- stats_summary_long %>%
  dplyr::filter(stringr::str_detect(quadrant, "TA\\+")) %>%
  dplyr::mutate(
    quadrant = dplyr::if_else(
      stringr::str_detect(quadrant, "not"),
      "Not targeted",
      "Targeted"
    )
  )

annotations_TA <- annotations %>%
  dplyr::filter(stringr::str_detect(quadrant, "TA\\+")) %>%
  dplyr::mutate(
    quadrant = dplyr::if_else(
      stringr::str_detect(quadrant, "not"),
      "Not targeted",
      "Targeted"
    )
  )

lines_TA <- lines %>%
  dplyr::filter(stringr::str_detect(quadrant, "TA\\+")) %>%
  dplyr::mutate(
    quadrant = dplyr::if_else(
      stringr::str_detect(quadrant, "not"),
      "Not targeted",
      "Targeted"
    )
  )

ticks_TA <- ticks %>%
  dplyr::filter(stringr::str_detect(quadrant, "TA\\+")) %>%
  dplyr::mutate(
    quadrant = dplyr::if_else(
      stringr::str_detect(quadrant, "not"),
      "Not targeted",
      "Targeted"
    )
  )

plot_02 <- stats_summary_long_TA %>%
  ggplot2::ggplot(ggplot2::aes(
    x = factor(group, x_order),
    y = n,
    fill = direction
  )) +
  ggplot2::geom_col(
    colour = "black",
    linewidth = 0.24,
    position = ggplot2::position_dodge(0.9)
  ) +
  ggplot2::scale_fill_manual(values = c("#e5e5e9", "#dc6465", "#5496ce")) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = lower, ymax = upper),
    position = ggplot2::position_dodge(0.9),
    width = 0.2,
    linewidth = 0.24
  ) +
  ggtext::geom_richtext(
    data = annotations_TA,
    ggplot2::aes(x = group, y = label_y, label = label),
    inherit.aes = FALSE,
    fill = NA,
    label.color = NA,
    size = 5 / ggplot2::.pt,
    hjust = 0.5
  ) +
  ggplot2::geom_segment(
    data = lines_TA,
    ggplot2::aes(x = x_start, xend = x_end, y = y, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.24,
    color = "black"
  ) +
  ggplot2::geom_segment(
    data = ticks_TA,
    ggplot2::aes(x = x, xend = x, y = y, yend = yend),
    inherit.aes = FALSE,
    linewidth = 0.24,
    color = "black"
  ) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(mult = c(0, 0.1))
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 0,
      hjust = 0.5,
      size = 9,
      colour = "black"
    ),
    legend.position = "top",
    text = ggplot2::element_text(size = 9, colour = "black"),
    axis.text = ggplot2::element_text(size = 9, colour = "black"),
    line = ggplot2::element_line(linewidth = 72.27 / 96 * 0.5),
    rect = ggplot2::element_rect(linewidth = 72.27 / 96 * 0.5),
    panel.border = ggplot2::element_rect(linewidth = 72.27 / 96 * 0.5),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.minor.x = ggplot2::element_blank(),
    panel.grid.minor.y = ggplot2::element_blank()
  ) +
  ggplot2::labs(
    y = "Number of plasmids carrying TA",
    x = "Type of targeting CRISPR-Cas"
  ) +
  ggplot2::coord_cartesian(clip = "off") +
  ggplot2::facet_wrap(~quadrant, scales = "free")

plot_02

plot_02 %>%
  ggplot2::ggsave(
    filename = "output/expected_vs_observed_02.pdf",
    width = 182.4,
    height = 80,
    units = "mm"
  )

plasmid_targeting_network_filt <- plasmid_targeting_network %>%
  dplyr::filter(
    source %in%
      plasmid_representatives_nover &
      target %in% plasmid_representatives_nover &
      source != target
  )

plasmid_targeting_network_filt_no_iv <- plasmid_targeting_network %>%
  dplyr::filter(
    source %in% plasmids_w_cas_no_iv_nover & source != target
  )

plasmid_targeting_network_filt_only_iv <- plasmid_targeting_network %>%
  dplyr::filter(
    source %in% plasmids_w_cas_only_iv_nover & source != target
  )

targeted_plasmids <- plasmid_targeting_network_filt %>%
  dplyr::distinct(target) %>%
  dplyr::pull()

targeted_plasmids_no_iv <- plasmid_targeting_network_filt_no_iv %>%
  dplyr::distinct(target) %>%
  dplyr::pull()

targeted_plasmids_only_iv <- plasmid_targeting_network_filt_only_iv %>%
  dplyr::distinct(target) %>%
  dplyr::pull()

tmp <- dplyr::bind_rows(
  dplyr::mutate(
    plasmid_targeting_network_filt,
    has_ta = dplyr::if_else(target %in% plasmids_w_ta_nover, TRUE, FALSE),
    category = "All"
  ),
  dplyr::mutate(
    plasmid_targeting_network_filt_no_iv,
    has_ta = dplyr::if_else(target %in% plasmids_w_ta_nover, TRUE, FALSE),
    category = "Type IV"
  ),
  dplyr::mutate(
    plasmid_targeting_network_filt_only_iv,
    has_ta = dplyr::if_else(target %in% plasmids_w_ta_nover, TRUE, FALSE),
    category = "Non-type IV"
  ),
) %>%
  dplyr::distinct(
    category,
    target,
    has_ta
  ) %>%
  dplyr::summarise(
    n = dplyr::n(),
    .by = c(category, has_ta)
  ) %>%
  dplyr::mutate(p = n / sum(n), .by = category) %>%
  dplyr::filter(has_ta == TRUE)

tmp


tmp %>%
  ggplot2::ggplot(ggplot2::aes(x = category, y = p)) +
  ggplot2::geom_col(colour = "black", fill = "#e5e5e9", linewidth = 0.24) +
  ggplot2::scale_y_continuous(
    limits = c(0, 1),
    expand = ggplot2::expansion(),
    labels = scales::label_percent()
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 45,
      hjust = 1,
      size = 5,
      colour = "black"
    ),
    legend.position = "top",
    text = ggplot2::element_text(size = 9, colour = "black"),
    axis.text = ggplot2::element_text(size = 9, colour = "black"),
    line = ggplot2::element_line(linewidth = 72.27 / 96 * 0.5),
    rect = ggplot2::element_rect(linewidth = 72.27 / 96 * 0.5),
    panel.border = ggplot2::element_rect(linewidth = 72.27 / 96 * 0.5),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.minor.x = ggplot2::element_blank(),
    panel.grid.minor.y = ggplot2::element_blank()
  ) +
  ggplot2::labs(
    y = "Proportion of plasmids carrying TA",
    x = "Type of targeting CRISPR-Cas"
  ) +
  ggplot2::coord_cartesian(clip = "off")

# TA ON CRISPR-BEARING PLASMIDS

ta_candidates_extended <- hits_best %>%
  dplyr::mutate(
    component = dplyr::if_else(
      stringr::str_detect(query, "^T"),
      "toxin",
      "antitoxin"
    )
  ) %>%
  dplyr::mutate(
    candidate = dplyr::case_when(
      "toxin" %in% component & "antitoxin" %in% component ~ "Complete TA",
      "toxin" %in% component & !"antitoxin" %in% component ~ "Orphan Toxin",
      !"toxin" %in% component & "antitoxin" %in% component ~ "Orphan Antitoxin"
    ),
    .by = c(seqid, family_cluster, family)
  )

hits_antitoxin <- ta_candidates_extended |>
  dplyr::filter(stringr::str_detect(query, "^AT")) |>
  dplyr::filter(candidate == "Complete TA" | candidate == "Orphan Antitoxin") |>
  dplyr::select(seqid, ta_type, family, candidate)

complete_ta <- hits_best |>
  dplyr::filter(stringr::str_detect(query, "^AT")) |>
  dplyr::distinct(subject, ta_type, family) |>
  dplyr::mutate(orphan = FALSE)

summary_antitoxins <- hits_antitoxin |>
  dplyr::distinct()

antitoxin_types_by_plasmid <- summary_antitoxins |>
  dplyr::rename(plasmid_id = seqid) |>
  dplyr::mutate(
    plasmid_id_nover = stringr::str_remove_all(plasmid_id, "\\..*")
  ) |>
  dplyr::mutate(
    at_type = list(ta_type),
    at_family = list(family),
    candidate = list(candidate),
    .by = plasmid_id
  ) |>
  dplyr::select(plasmid_id, plasmid_id_nover, at_type, at_family, candidate) |>
  dplyr::distinct()

ta_types_by_plasmid <- summary_antitoxins |>
  dplyr::filter(candidate == "Complete TA") |>
  dplyr::rename(plasmid_id = seqid) |>
  dplyr::mutate(
    plasmid_id_nover = stringr::str_remove_all(plasmid_id, "\\..*")
  ) |>
  dplyr::mutate(
    ta_type = list(ta_type),
    ta_family = list(family),
    candidate = list(candidate),
    .by = plasmid_id
  ) |>
  dplyr::select(plasmid_id, plasmid_id_nover, ta_type, ta_family) |>
  dplyr::distinct()

tmp <- plasmid_targeting_network_filt_no_iv |>
  dplyr::left_join(
    ta_types_by_plasmid,
    by = dplyr::join_by(target == plasmid_id_nover)
  ) |>
  # dplyr::filter(!is.na(plasmid_id)) |>
  dplyr::select(!plasmid_id) |>
  dplyr::rename(target_ta_types = ta_type, target_ta_families = ta_family) |>
  dplyr::left_join(
    ta_types_by_plasmid,
    by = dplyr::join_by(source == plasmid_id_nover)
  ) |>
  # dplyr::filter(!is.na(plasmid_id)) |>
  dplyr::select(!plasmid_id) |>
  dplyr::rename(source_ta_types = ta_type, source_ta_families = ta_family) |>
  dplyr::distinct()

tmp2 <- tmp |>
  dplyr::filter(lengths(target_ta_types) > 0) |>
  dplyr::mutate(
    overlapping_ta = purrr::map2(
      target_ta_families,
      source_ta_families,
      intersect
    )
  ) |>
  dplyr::mutate(
    compatible_overlap = purrr::map2_lgl(
      target_ta_families,
      source_ta_families,
      ~ length(setdiff(.x, .y)) == 0
    )
  ) |>
  dplyr::mutate(partial_overlap = lengths(overlapping_ta) > 0)

non_type_IV_ta_compatible_overlap_summary <- tmp2 |>
  dplyr::summarise(n = dplyr::n(), .by = compatible_overlap) |>
  dplyr::mutate(total = sum(n)) |>
  dplyr::mutate(p = n / total * 100)

non_type_IV_ta_compatible_overlap_summary

non_type_IV_ta_partial_overlap_summary <- tmp2 |>
  dplyr::summarise(n = dplyr::n(), .by = partial_overlap) |>
  dplyr::mutate(total = sum(n)) |>
  dplyr::mutate(p = n / total * 100)

non_type_IV_ta_partial_overlap_summary

tmp3 <- plasmid_targeting_network_filt_no_iv |>
  dplyr::left_join(
    ta_types_by_plasmid,
    by = dplyr::join_by(target == plasmid_id_nover)
  ) |>
  # dplyr::filter(!is.na(plasmid_id)) |>
  dplyr::select(!plasmid_id) |>
  dplyr::rename(target_ta_types = ta_type, target_ta_families = ta_family) |>
  dplyr::left_join(
    antitoxin_types_by_plasmid,
    by = dplyr::join_by(source == plasmid_id_nover)
  ) |>
  # dplyr::filter(!is.na(plasmid_id)) |>
  dplyr::select(!plasmid_id) |>
  dplyr::rename(source_at_types = at_type, source_at_families = at_family) |>
  dplyr::distinct()

tmp4 <- tmp3 |>
  dplyr::filter(lengths(target_ta_types) > 0) |>
  dplyr::mutate(
    overlapping_at = purrr::map2(
      target_ta_families,
      source_at_families,
      intersect
    )
  ) |>
  dplyr::mutate(
    compatible_overlap = purrr::map2_lgl(
      target_ta_families,
      source_at_families,
      ~ length(setdiff(.x, .y)) == 0
    )
  ) |>
  dplyr::mutate(partial_overlap = lengths(overlapping_at) > 0)

non_type_IV_at_compatible_overlap_summary <- tmp4 |>
  dplyr::summarise(compatible_overlap = dplyr::n(), .by = compatible_overlap) |>
  dplyr::mutate(total = sum(compatible_overlap)) |>
  dplyr::mutate(p = compatible_overlap / total * 100)

non_type_IV_at_compatible_overlap_summary

### Type IV

tmp5 <- plasmid_targeting_network_filt_only_iv |>
  dplyr::left_join(
    ta_types_by_plasmid,
    by = dplyr::join_by(target == plasmid_id_nover)
  ) |>
  # dplyr::filter(!is.na(plasmid_id)) |>
  dplyr::select(!plasmid_id) |>
  dplyr::rename(target_ta_types = ta_type, target_ta_families = ta_family) |>
  dplyr::left_join(
    ta_types_by_plasmid,
    by = dplyr::join_by(source == plasmid_id_nover)
  ) |>
  # dplyr::filter(!is.na(plasmid_id)) |>
  dplyr::select(!plasmid_id) |>
  dplyr::rename(source_ta_types = ta_type, source_ta_families = ta_family) |>
  dplyr::distinct() |>
  dplyr::filter(lengths(target_ta_types) > 0)

tmp6 <- tmp5 |>
  dplyr::mutate(
    overlapping_ta = purrr::map2(
      target_ta_families,
      source_ta_families,
      intersect
    )
  ) |>
  dplyr::mutate(
    compatible_overlap = purrr::map2_lgl(
      target_ta_families,
      source_ta_families,
      ~ length(setdiff(.x, .y)) == 0
    )
  ) |>
  dplyr::mutate(partial_overlap = lengths(overlapping_ta) > 0)

only_type_IV_ta_compatible_overlap_summary <- tmp6 |>
  dplyr::summarise(n = dplyr::n(), .by = compatible_overlap) |>
  dplyr::mutate(total = sum(n)) |>
  dplyr::mutate(p = n / total * 100)

only_type_IV_ta_compatible_overlap_summary

only_type_IV_ta_partial_overlap_summary <- tmp6 |>
  dplyr::summarise(n = dplyr::n(), .by = partial_overlap) |>
  dplyr::mutate(total = sum(n)) |>
  dplyr::mutate(p = n / total * 100)

only_type_IV_ta_partial_overlap_summary

tmp7 <- plasmid_targeting_network_filt_only_iv |>
  dplyr::left_join(
    ta_types_by_plasmid,
    by = dplyr::join_by(target == plasmid_id_nover)
  ) |>
  # dplyr::filter(!is.na(plasmid_id)) |>
  dplyr::select(!plasmid_id) |>
  dplyr::rename(target_ta_types = ta_type, target_ta_families = ta_family) |>
  dplyr::left_join(
    antitoxin_types_by_plasmid,
    by = dplyr::join_by(source == plasmid_id_nover)
  ) |>
  # dplyr::filter(!is.na(plasmid_id)) |>
  dplyr::select(!plasmid_id) |>
  dplyr::rename(source_at_types = at_type, source_at_families = at_family) |>
  dplyr::distinct()

tmp8 <- tmp7 |>
  dplyr::filter(lengths(target_ta_types) > 0) |>
  dplyr::mutate(
    overlapping_at = purrr::map2(
      target_ta_families,
      source_at_families,
      intersect
    )
  ) |>
  dplyr::mutate(
    compatible_overlap = purrr::map2_lgl(
      target_ta_families,
      source_at_families,
      ~ length(setdiff(.x, .y)) == 0
    )
  ) |>
  dplyr::mutate(partial_overlap = lengths(overlapping_at) > 0)

only_type_IV_at_compatible_overlap_summary <- tmp8 |>
  dplyr::summarise(compatible_overlap = dplyr::n(), .by = compatible_overlap) |>
  dplyr::mutate(total = sum(compatible_overlap)) |>
  dplyr::mutate(p = compatible_overlap / total * 100)

only_type_IV_at_compatible_overlap_summary


# STATS

non_type_IV_ta_compatibility <- tmp3 |>
  dplyr::mutate(
    overlapping_at = purrr::map2(
      target_ta_families,
      source_at_families,
      intersect
    )
  ) |>
  dplyr::mutate(
    compatible_overlap = purrr::map2_lgl(
      target_ta_families,
      source_at_families,
      ~ length(setdiff(.x, .y)) == 0
    )
  ) |>
  dplyr::mutate(
    target_has_ta = dplyr::case_when(
      lengths(target_ta_types) > 0 ~ TRUE,
      .default = FALSE
    ),
    targeting_compatible = dplyr::case_when(
      compatible_overlap == TRUE ~ TRUE,
      .default = FALSE
    )
  ) |>
  dplyr::select(source, target, target_has_ta, targeting_compatible)

non_type_IV_ta_compatibility

only_type_IV_ta_compatibility <- tmp7 |>
  dplyr::mutate(
    overlapping_at = purrr::map2(
      target_ta_families,
      source_at_families,
      intersect
    )
  ) |>
  dplyr::mutate(
    compatible_overlap = purrr::map2_lgl(
      target_ta_families,
      source_at_families,
      ~ length(setdiff(.x, .y)) == 0
    )
  ) |>
  dplyr::mutate(
    target_has_ta = dplyr::case_when(
      lengths(target_ta_types) > 0 ~ TRUE,
      .default = FALSE
    ),
    targeting_compatible = dplyr::case_when(
      compatible_overlap == TRUE ~ TRUE,
      .default = FALSE
    )
  ) |>
  dplyr::select(source, target, target_has_ta, targeting_compatible)

only_type_IV_ta_compatibility

all_ta_compatibility <-
  dplyr::bind_rows(
    non_type_IV_ta_compatibility |> dplyr::mutate(type = "nucleolytic"),
    only_type_IV_ta_compatibility |> dplyr::mutate(type = "non-nucleolytic")
  )

all_ta_compatibility

all_ta_compatibility_where_target_has_ta <- all_ta_compatibility |>
  dplyr::filter(target_has_ta)

# PLOTTING #1

contingency_table_all <- table(
  all_ta_compatibility_where_target_has_ta$type,
  all_ta_compatibility_where_target_has_ta$targeting_compatible
)

ft_all <- fisher.test(contingency_table_all)
ft_all

N <- sum(contingency_table_all)
rs <- rowSums(contingency_table_all) # row totals
cs <- colSums(contingency_table_all) # col totals

# full expected matrix under independence
E <- outer(rs, cs) / N

# exact 95% null interval under hypergeometric distribution
n <- as.numeric(rs) # row totals (vector)
K <- as.numeric(cs["TRUE"]) # col total for "compatible=TRUE"
lower_null <- qhyper(0.025, K, N - K, n)
upper_null <- qhyper(0.975, K, N - K, n)

# Haberman (standardized) residuals
residuals_matrix_all <- (contingency_table_all - E) /
  sqrt(E * (1 - rs / N) %o% (1 - cs / N))
residuals_matrix_all

# Tidy summary just for "compatible = TRUE"
e_vs_o_all <- tibble::tibble(
  group = rownames(contingency_table_all),
  quadrant = "Compatible antitoxin",
  observed_count = as.numeric(contingency_table_all[, "TRUE"]),
  expected_count = as.numeric(E[, "TRUE"]),
  lower = lower_null,
  upper = upper_null,
  residual = as.numeric(residuals_matrix_all[, "TRUE"]),
  p_value = ft_all$p.value,
  direction = ifelse(
    observed_count > expected_count,
    "Observed (+)",
    "Observed (-)"
  )
)

e_vs_o_all

x_order <- c("non-nucleolytic", "nucleolytic")

stats_summary_long <- e_vs_o_all %>%
  tidyr::pivot_longer(
    c(observed_count, expected_count),
    names_to = "category",
    values_to = "n"
  ) %>%
  dplyr::mutate(
    category = dplyr::if_else(
      category == "observed_count",
      "Observed",
      "Expected"
    ),
    direction = dplyr::if_else(category == "Observed", direction, "Expected"),
    lower = ifelse(category == "Expected", lower, NA),
    upper = ifelse(category == "Expected", upper, NA)
  )

format_p <- function(p) {
  ifelse(
    p > 0.001,
    formatC(p, format = "f", digits = 4),
    gsub(
      "e\\-0*(\\d+)",
      " × 10<sup>−\\1</sup>",
      formatC(p, format = "e", digits = 2)
    )
  )
}

annotations <- stats_summary_long %>%
  dplyr::filter(category == "Observed") %>%
  dplyr::mutate(
    label_y = n * 1.2,
    formatted_p = format_p(p_value),
    label = paste0(
      "Residual = ",
      round(residual, 2),
      "<br><i>p</i> = ",
      formatted_p
    )
  ) %>%
  dplyr::select(group, quadrant, label, label_y)

annotations

plot_fisher <- stats_summary_long %>%
  ggplot2::ggplot(
    ggplot2::aes(
      x = factor(group, x_order),
      y = n,
      fill = direction
    )
  ) +
  ggplot2::geom_col(
    colour = "black",
    linewidth = 0.24,
    position = ggplot2::position_dodge(0.9)
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      "Expected" = "#e5e5e9",
      "Observed (+)" = "#ad4c19",
      "Observed (-)" = "#ad4c19"
    )
  ) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = lower, ymax = upper),
    position = ggplot2::position_dodge(0.9),
    width = 0.2,
    linewidth = 0.24
  ) +
  ggtext::geom_richtext(
    data = annotations,
    ggplot2::aes(x = group, y = label_y, label = label),
    inherit.aes = FALSE,
    fill = NA,
    label.color = NA,
    size = 5 / ggplot2::.pt,
    hjust = 0.5
  ) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(mult = c(0, 0.1))
  ) +
  ggplot2::labs(
    y = "Count with compatible antitoxin (TA+ targets)",
    x = "Type of targeting CRISPR-Cas"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 0,
      hjust = 0.5,
      size = 9,
      colour = "black"
    ),
    legend.position = "none",
    text = ggplot2::element_text(size = 9, colour = "black"),
    axis.text = ggplot2::element_text(size = 9, colour = "black"),
    line = ggplot2::element_line(linewidth = 72.27 / 96 * 0.5),
    rect = ggplot2::element_rect(linewidth = 72.27 / 96 * 0.5),
    panel.border = ggplot2::element_rect(linewidth = 72.27 / 96 * 0.5),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.minor.x = ggplot2::element_blank(),
    panel.grid.minor.y = ggplot2::element_blank(),
  )

plot_fisher

ggplot2::ggsave(
  plot_fisher,
  filename = "output/ta_compatibility.pdf",
  width = 91.2,
  height = 60,
  units = "mm"
)

sup <- dplyr::bind_rows(
  tmp3 |> dplyr::mutate(`Source CRISPR-Cas Type` = "Non-type IV"),
  tmp7 |> dplyr::mutate(`Source CRISPR-Cas Type` = "Type IV")
) |>
  dplyr::mutate(
    overlapping_at = purrr::map2(
      target_ta_families,
      source_at_families,
      intersect
    )
  ) |>
  dplyr::mutate(
    compatible = purrr::map2_lgl(
      target_ta_families,
      source_at_families,
      ~ length(setdiff(.x, .y)) == 0
    )
  ) |>
  dplyr::rowwise() |>
  dplyr::mutate(
    Compatible = dplyr::if_else(
      length(target_ta_types) > 0,
      as.character(compatible),
      ""
    )
  ) |>
  dplyr::select(!c(overlapping_at, compatible)) |>
  dplyr::rowwise() |>
  dplyr::mutate(dplyr::across(
    dplyr::everything(),
    ~ paste0(., collapse = ", ")
  )) |>
  dplyr::rename(
    Source = source,
    Target = target,
    `Target TA Types` = target_ta_types,
    `Target TA Families` = target_ta_families,
    `Source AT Types` = source_at_types,
    `Source AT Families` = source_at_families,
    `Source AT Categories` = candidate
  ) |>
  dplyr::select(
    `Source`,
    `Target`,
    `Source CRISPR-Cas Type`,
    dplyr::everything()
  ) |>
  dplyr::arrange(
    Source,
    Target,
    `Source CRISPR-Cas Type`
  )

sup |>
  writexl::write_xlsx("output/ta_compatibility.xlsx")

stats_out <- stats_summary_long |>
  dplyr::select(
    group,
    category,
    n,
    lower,
    upper,
    residual,
    p_value
  ) |>
  dplyr::rename(count = n) |>
  dplyr::mutate(
    group = dplyr::case_when(
      group == "non-nucleolytic" ~ "Type IV",
      .default = "Non-type IV"
    )
  )

stats_out

stats_out |> 
  writexl::write_xlsx("output/ta_compatibility_stats.xlsx")
