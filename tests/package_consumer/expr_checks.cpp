int run_expr_arithmetic_consumer_checks();
int run_expr_complex_consumer_checks();
int run_expr_sets_consumer_checks();
int run_expr_contracts_consumer_checks();

int run_expr_consumer_checks() {
    if (const int status = run_expr_arithmetic_consumer_checks(); status != 0) {
        return status;
    }
    if (const int status = run_expr_complex_consumer_checks(); status != 0) {
        return status;
    }
    if (const int status = run_expr_sets_consumer_checks(); status != 0) {
        return status;
    }
    if (const int status = run_expr_contracts_consumer_checks(); status != 0) {
        return status;
    }
    return 0;
}
