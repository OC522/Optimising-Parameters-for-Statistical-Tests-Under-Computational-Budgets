# ============================================================
# Run Chapter 5 experiments using runtime-fixed algorithm3.1
# ============================================================
chapter5_file <- if (file.exists("chapter5_runtime_fixed.R")) {
  "chapter5_runtime_fixed.R"
} else {
  "~/Downloads/chapter5_runtime_fixed.R"
}
source(chapter5_file)

res_ch5 <- run_ch5_all_runtime_fixed()

cat("\nBudget diagnostics: main\n")
print(budget_diagnostics(res_ch5$main))
cat("\nBudget diagnostics: ablation\n")
print(budget_diagnostics(res_ch5$ablation))
cat("\nBudget diagnostics: lambda\n")
print(budget_diagnostics(res_ch5$lambda))
