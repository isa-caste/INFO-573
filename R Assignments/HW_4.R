## Task 1
# load in the file 
data = readLines("C:/Users/icwin/OneDrive - Indiana University/[B-573] Prog. Sci-Informatics/Python/chr1_GL383518v1_alt.fa")[-1]
# Print the 10th letter of this sequence
substring <- paste(data, collapse = "")
cat(substr(substring, 10, 10))
# Print the 758th letter of this sequence
cat(substr(substring, 758, 758))

## Task 2
# create function to create the reverse compliment string
reverse_complement <- function(dna_seq) {
  complement <- chartr("ATCG", "TAGC", dna_seq)
  reverse_complement_seq <- rev(strsplit(complement, NULL)[[1]])
  return(paste(reverse_complement_seq, collapse = ""))
}
rev_comp_sequence <- reverse_complement(substring)
# print the 79th letter of this sequence.
cat("79th letter:", substr(rev_comp_sequence, 79, 79), "\n")
# print the 500th through the 800th letters of this sequence.
cat("500th to 800th letters:\n", substr(rev_comp_sequence, 500, 800), "\n")

## Task 3
# create function to count each specific nucleotide
count_nucleotides_per_kilobase <- function(dna_seq) {
  num_kilobases <- ceiling(nchar(dna_seq) / 1000)
  # create empty list for counting the bases
  kilobase_counts <- vector("list", num_kilobases)
  for (i in 1:num_kilobases) {
    start <- (i - 1) * 1000 + 1
    end <- min(i * 1000, nchar(dna_seq))
    segment <- substr(dna_seq, start, end)
    # start counting each base pair
    counts <- list(
      A = sum(strsplit(segment, NULL)[[1]] == "A"),
      T = sum(strsplit(segment, NULL)[[1]] == "T"),
      C = sum(strsplit(segment, NULL)[[1]] == "C"),
      G = sum(strsplit(segment, NULL)[[1]] == "G")
    )
    kilobase_counts[[i]] <- counts
  }
  return(kilobase_counts)
}
print(count_nucleotides_per_kilobase(substring))

## Task 4
# create function to count nucleotides per kilobase
count_nucleotides_per_kilobase <- function(dna_seq) {
  num_kilobases <- ceiling(nchar(dna_seq) / 1000)
  kilobase_counts <- vector("list", num_kilobases)
  # loop through to find kilobases
  for (i in 1:num_kilobases) {
    start <- (i - 1) * 1000 + 1
    end <- min(i * 1000, nchar(dna_seq))
    segment <- substr(dna_seq, start, end)
    # count nucleotides
    counts <- list(
      A = sum(strsplit(segment, NULL)[[1]] == "A"),
      T = sum(strsplit(segment, NULL)[[1]] == "T"),
      C = sum(strsplit(segment, NULL)[[1]] == "C"),
      G = sum(strsplit(segment, NULL)[[1]] == "G")
    )
    kilobase_counts[[i]] <- counts
  }
  return(kilobase_counts)
}
print(count_nucleotides_per_kilobase(substring))


