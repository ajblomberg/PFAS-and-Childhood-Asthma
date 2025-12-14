# 00_run_all_public.R
# Purpose: Run the full PFAS-asthma analysis in a fixed, documented order
# Usage: source ("00_run_all_public.R") 

# 1. Basic run configuration -------------------------------------------------

set.seed(1516)

code_dir <- "code_public"
output_dir <- "outputs_public"

if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

log_path <- file.path(output_dir, "run_all.log")


# 2. Start Log ------------------------------------------------

# Open log connection 
log_con <- file(log_path, open = "wt")

# Route both stdout + messages into the log 
sink(log_con, type = c("output", "message"), split = T)

# Close sinks and connections, even if an error occurs mid-run:
on.exit({
  sink(type = "message")
  sink(type = "output")
  close(log_con)
}, add = TRUE)

# Write run header into log
cat(sprintf("=== PFAS-asthma run started: %s === \n", Sys.time()))
cat("Project root: ", normalizePath("."), "\n")
cat(sprintf("Log file: %s\n", normalizePath(log_path)), "\n")
cat("R version ", R.version.string, "\n\n")


# 4. Record environment and session info ----------------------------------

sess_file <- file.path(output_dir, sprintf("sessionInfo_%s.txt", 
                                           format(Sys.Date(), "%Y-%m-%d")))
writeLines(c(capture.output(sessionInfo()), ""), con = sess_file)
cat("Saved session info to: ", sess_file, "\n\n")


# 5. Analysis run ------------------------------------------------------------

scripts <- list.files(code_dir, pattern = "^\\d+[-_].*\\.R$", full.names = TRUE)
if(length(scripts) == 0L) stop("No numbered scripts found in ./code")

# Order by the numeric prefix, not lexicographically
num_from_basename <- function(p) as.integer(sub("^([0-9]+).*", "\\1", basename(p)))
ord <- order(num_from_basename(scripts), basename(scripts))
scripts <- scripts[ord]

cat("Execution plan: \n")
for (p in scripts) cat(" - ", basename(p), "\n")
cat("\n")


# 6. Source each script with timing info -------------------------------------
run_times <- data.frame(script = basename(scripts), seconds = NA_real_)

for(i in seq_along(scripts)) {
  scr <- scripts[i]
  cat(sprintf("\n>>> [%02d/%02d] Running %s\n", i, length(scripts), basename(scr)))
  cat("Start time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  
  t <- system.time({
    # Use local() so objects created inside scripts don't pollute the global env by accident
    tryCatch(
      source(scr, echo = TRUE, max.deparse.length = Inf, local = TRUE),
      error = function(e) {
        cat("\nError in ", basename(scr), ":\n", conditionMessage(e), "\n", sep = "")
        stop(e)
      }
    )
    
  })
  
  run_times$seconds[i] <- unname(t["elapsed"])
  cat("Elapsed (s): ", sprintf("%.2f", run_times$seconds[i]), "\n")
}


# 7. Write simple run manifest and close -------------------------------------

# Save the timing table for transparency 
timing_csv <- file.path(output_dir, "run_times.csv")
utils::write.csv(run_times, timing_csv, row.names = FALSE)
cat("save run times to: ", timing_csv, "\n")

cat(sprintf("\n== Run completed successfully: %s ===\n", Sys.time()))



