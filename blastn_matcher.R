### blastn OUTPUT TABLE PARSER AND SUCCESS RATE CALCULATOR
### 24.05.2019

library(tidyverse)
library(stringr)

# DATA LOADING ------------------------------------------------------------


# 1. OPENING THE BLASTN RESULTS TABLE (outfmt = 6); ONLY COLUMNS 1,2 3, AND 12 ARE SELECTED. SUBSITUTE "barcode" by "rbcLa", "matK" or 
# "trnH-psbA".
data.df <- read.csv("docs/barcode.blastn.SWFRG", sep = "\t", header = F) %>% select(V1, V2, V3, V12)

# 2. SETTING THE ACTUAL NAME FOR THE COLUMNS
colnames(data.df) <- c("QUERY", "REF", "PID", "SCORE")

# 3. COLUMN "QUERY" IS SEPARATED INTO THREE, IN ORDER TO HAVE THE SEQUENCE ID., THE QUERY FAMILY,
#   AND THE QUERY SPECIES IN SEPARATE COLUMNS. 
data.df <- data.df %>% separate(
  QUERY, sep = '\\|', into = c("QUERY.SEQID", "QUERY.FAMILY", "QUERY.SPECIES") )

# 4. SIMILARLY, COLUMN "REF" IS ALSO SEPARATED INTO THREE SEPARATE COLUMNS. 
#     SOME WARNINGS MAY ALSO APPEAR DUE TO THE separate() FUNCTION.
data.df <- data.df %>% separate(
  REF, sep = '\\|', into = c("REF.SEQID", "REF.FAMILY", "REF.SPECIES")) 

# 5. FOR OUR "LEAVE-ONE-OUT" BLAST APPROACH, WE BLASTED THE SAME FASTA FILE THAT WAS USED TO PRODUCE 
#   THE BLAST DATABASE. THIS MEANS THAT EACH FASTA ENTRY PRODUCED A HIT AGAINST ITSELF.
# 
#   IN THIS STEP WE REMOVE THOSE "SELF-HITS" BY FILTERING OUT THE HITS FOR WHICH THE QUERY SEQUENCE ID
#   IS THE SAME AS THE REFERENCE SEQUENCE ID.
data.df <- data.df %>% filter(QUERY.SEQID != REF.SEQID)

# 6. THERE ARE SOME TYPOS IN THE SPECIES NAMES FROM THE BOLD DATABASE. THESE ARE SOME WE FOUND:
data.df$QUERY.SPECIES <- gsub('Arrhenaterum_elatius', 'Arrhenatherum_elatius', data.df$QUERY.SPECIES)
data.df$QUERY.SPECIES <- gsub('Festuca_pratense', 'Festuca_pratensis', data.df$QUERY.SPECIES)
data.df$REF.SPECIES <- gsub('Arrhenaterum_elatius', 'Arrhenatherum_elatius', data.df$REF.SPECIES)
data.df$REF.SPECIES <- gsub('Festuca_pratense', 'Festuca_pratensis', data.df$REF.SPECIES)
data.df$QUERY.SPECIES <- gsub('Triflorum', 'Trifolium', data.df$QUERY.SPECIES)
data.df$REF.SPECIES <- gsub('Triflorum', 'Trifolium', data.df$REF.SPECIES)

# MATCH CALULATOR: GENERAL ------------------------------------------------
### A. GENERAL RESULTS
data.df %>% 
  # 7. FILTERING IN ONLY THE SAMPLES FROM OUR PROJECT: "SWFRG"
  filter("SWFRG" == str_sub(QUERY.SEQID, 1,5)) %>%
  
  # 8. FILTERING IN ONLY THE BEST SCORING HITS FOR EACH QUERY SEQUENCE
  group_by(QUERY.SEQID) %>%
  filter(SCORE == max(SCORE, na.rm = T)) %>%
  
  # 9. CALLING MATCHES WITH THE FOLLOWING CRITERIA:
  #   -IF THERE IS MORE THAN ONE HIGH-SCORING REF.SPECIES, IT'S A NON-MATCH
  #   -IF THE QUERY SPECIES HITS A UNIQUE REF.SPECIES, IT'S A MATCH
  mutate(
    MATCH = ifelse(length(unique(REF.SPECIES)) > 1,F,
                   ifelse(QUERY.SPECIES == unique(REF.SPECIES), T,F))) %>% 
  
  # 10. THIS STEPS SUMMARISES THE NUMBER OF HITS, INCLUDING MATCHES AND NON-MATCHES, 
  #     FROM EACH QUERY SEQUENCE IN ONE ROW
  group_by(QUERY.SEQID, QUERY.SPECIES, MATCH, QUERY.FAMILY) %>% summarise(n()) %>% 
  
  # 11. FINALLY, THIS STEP SUMMARISES THE NUMBER OF MATCHES OR NON-MATCHES, AS WELL
  #     AS THEIR FREQUENCIES
  group_by(MATCH) %>%
  summarise(n = n()) %>% 
  mutate(freq = n/sum(n))

### B. RESULTS BY PLANT FAMILY
data.df %>% 
  # 12. FILTERING IN ONLY THE SAMPLES FROM OUR PROJECT: "SWFRG"
  filter("SWFRG" == str_sub(QUERY.SEQID, 1,5)) %>%
  
  # 13. FILTERING IN ONLY THE BEST SCORING HITS FOR EACH QUERY SEQUENCE
   group_by(QUERY.SEQID) %>%
  filter(SCORE == max(SCORE, na.rm = T)) %>% 
  
  # 14. CALLING MATCHES WITH THE FOLLOWING CRITERIA:
  #   -IF THERE IS MORE THAN ONE HIGH-SCORING REF.SPECIES, IT'S A NON-MATCH
  #   -IF THE QUERY SPECIES HITS A UNIQUE REF.SPECIES, IT'S A MATCH
  mutate(
    MATCH = ifelse(length(unique(REF.SPECIES)) > 1,F,
                   ifelse(QUERY.SPECIES == unique(REF.SPECIES), T,F))) %>% 
  
  # 15. THIS STEPS SUMMARISES THE NUMBER OF HITS, INCLUDING MATCHES AND NON-MATCHES, 
  #     FROM EACH QUERY SEQUENCE IN ONE ROW
  group_by(QUERY.SEQID, QUERY.SPECIES, MATCH, QUERY.FAMILY) %>% summarise(n()) %>% 
  
  # 16. FINALLY, THIS STEP SUMMARISES THE NUMBER OF MATCHES OR NON-MATCHES, AS WELL
  #     AS THEIR FREQUENCIES
  group_by(MATCH, QUERY.FAMILY) %>%
  summarise(n = n()) %>% group_by(QUERY.FAMILY) %>%
  mutate(freq = n/sum(n))


# GENUS-LEVEL ANALYSIS ----------------------------------------------------

data.df2 <- data.df %>% 
  separate(QUERY.SPECIES, sep = '_', into = c('QUERY.GENUS', 'QUERY.SPEP')) %>%
  separate(REF.SPECIES, sep = '_', into = c('REF.GENUS', 'REF.SPEP'))

data.df2 %>% 
  # 17. FILTERING IN ONLY THE SAMPLES FROM OUR PROJECT: "SWFRG"
  filter("SWFRG" == str_sub(QUERY.SEQID, 1,5)) %>%
  
  # 18. FILTERING IN ONLY THE BEST SCORING HITS FOR EACH QUERY SEQUENCE
  group_by(QUERY.SEQID) %>%
  filter(SCORE == max(SCORE, na.rm = T)) %>%
  
  # 19. CALLING MATCHES WITH THE FOLLOWING CRITERIA:
  #   -IF THERE IS MORE THAN ONE HIGH-SCORING REF.SPECIES, IT'S A NON-MATCH
  #   -IF THE QUERY SPECIES HITS A UNIQUE REF.SPECIES, IT'S A MATCH
  mutate(
    MATCH = ifelse(length(unique(REF.GENUS)) > 1,F,
                   ifelse(QUERY.GENUS == unique(REF.GENUS), T,F))) %>% 
  
  # 20. THIS STEPS SUMMARISES THE NUMBER OF HITS, INCLUDING MATCHES AND NON-MATCHES, 
  #     FROM EACH QUERY SEQUENCE IN ONE ROW
  group_by(QUERY.SEQID, QUERY.GENUS, MATCH, QUERY.FAMILY) %>% summarise(n()) %>% 
  
  # 21. FINALLY, THIS STEP SUMMARISES THE NUMBER OF MATCHES OR NON-MATCHES, AS WELL
  #     AS THEIR FREQUENCIES
  group_by(MATCH) %>%
  summarise(n = n()) %>% 
  mutate(freq = n/sum(n))

### B. RESULTS BY PLANT FAMILY
data.df2 %>% 
  # 22. FILTERING IN ONLY THE SAMPLES FROM OUR PROJECT: "SWFRG"
  filter("SWFRG" == str_sub(QUERY.SEQID, 1,5)) %>%
  
  # 23. FILTERING IN ONLY THE BEST SCORING HITS FOR EACH QUERY SEQUENCE
  group_by(QUERY.SEQID) %>%
  filter(SCORE == max(SCORE, na.rm = T)) %>%
  
  # 24. CALLING MATCHES WITH THE FOLLOWING CRITERIA:
  #   -IF THERE IS MORE THAN ONE HIGH-SCORING REF.SPECIES, IT'S A NON-MATCH
  #   -IF THE QUERY SPECIES HITS A UNIQUE REF.SPECIES, IT'S A MATCH
  mutate(
    MATCH = ifelse(length(unique(REF.GENUS)) > 1,F,
                   ifelse(QUERY.GENUS == unique(REF.GENUS), T,F))) %>% 
  
  # 25. THIS STEPS SUMMARISES THE NUMBER OF HITS, INCLUDING MATCHES AND NON-MATCHES, 
  #     FROM EACH QUERY SEQUENCE IN ONE ROW
  group_by(QUERY.SEQID, QUERY.GENUS, MATCH, QUERY.FAMILY) %>% summarise(n()) %>% 
  
  # 26. FINALLY, THIS STEP SUMMARISES THE NUMBER OF MATCHES OR NON-MATCHES, AS WELL
  #     AS THEIR FREQUENCIES
  group_by(MATCH, QUERY.FAMILY) %>%
  summarise(n = n()) %>% group_by(QUERY.FAMILY) %>%
  mutate(freq = n/sum(n))
