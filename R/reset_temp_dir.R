

# Function to reset temporary directory for each replication
reset_temp_dir <- function(tempdir_path) {
  # If the directory exists, delete it
  if (dir.exists(tempdir_path)) {
    unlink(tempdir_path, recursive = TRUE)  # Delete directory and its contents
  }

  # Recreate the temporary directory
  dir.create(tempdir_path, recursive = TRUE)

  # Set the directory for temporary files in R and bigstatsr
  Sys.setenv(TMPDIR = tempdir_path)
  options(bigstatsr.temporary_directory = tempdir_path)
}
