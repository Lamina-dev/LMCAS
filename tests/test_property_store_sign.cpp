
#include "test_common.hpp"
#include "property_store.hpp"
#include "symbolic.hpp"

using namespace LMCAS;

TEST(PropertyStoreSign, PositiveImpliesNonnegativeAndNonzero) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_sign("x", Sign::Positive).has_value())) << "sign declaration succeeds";

    EXPECT_TRUE((store.has_sign("x", Sign::Positive))) << "x has Positive";
    EXPECT_TRUE((store.has_sign("x", Sign::NonNegative))) << "x has NonNegative (implied)";
    EXPECT_TRUE((store.has_sign("x", Sign::NonZero))) << "x has NonZero (implied)";
    EXPECT_FALSE((store.has_sign("x", Sign::Negative))) << "x does not have Negative";
    EXPECT_FALSE((store.has_sign("x", Sign::NonPositive))) << "x does not have NonPositive";
    EXPECT_FALSE((store.has_sign("x", Sign::Zero))) << "x does not have Zero";
}

TEST(PropertyStoreSign, NegativeImpliesNonpositiveAndNonzero) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_sign("y", Sign::Negative).has_value())) << "sign declaration succeeds";

    EXPECT_TRUE((store.has_sign("y", Sign::Negative))) << "y has Negative";
    EXPECT_TRUE((store.has_sign("y", Sign::NonPositive))) << "y has NonPositive (implied)";
    EXPECT_TRUE((store.has_sign("y", Sign::NonZero))) << "y has NonZero (implied)";
    EXPECT_FALSE((store.has_sign("y", Sign::Positive))) << "y does not have Positive";
    EXPECT_FALSE((store.has_sign("y", Sign::NonNegative))) << "y does not have NonNegative";
    EXPECT_FALSE((store.has_sign("y", Sign::Zero))) << "y does not have Zero";
}

TEST(PropertyStoreSign, ZeroImpliesNonnegativeNonpositiveAndInteger) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_sign("z", Sign::Zero).has_value())) << "sign declaration succeeds";

    EXPECT_TRUE((store.has_sign("z", Sign::Zero))) << "z has Zero";
    EXPECT_TRUE((store.has_sign("z", Sign::NonNegative))) << "z has NonNegative (implied)";
    EXPECT_TRUE((store.has_sign("z", Sign::NonPositive))) << "z has NonPositive (implied)";
    EXPECT_FALSE((store.has_sign("z", Sign::Positive))) << "z does not have Positive";
    EXPECT_FALSE((store.has_sign("z", Sign::Negative))) << "z does not have Negative";
    EXPECT_FALSE((store.has_sign("z", Sign::NonZero))) << "z does not have NonZero";

    // Zero also implies Integer domain
    EXPECT_TRUE((store.has_domain("z", Domain::Integer))) << "z has Integer domain (implied by Zero)";
}

TEST(PropertyStoreSign, NonnegativeNoExtraImplications) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_sign("a", Sign::NonNegative).has_value())) << "sign declaration succeeds";

    EXPECT_TRUE((store.has_sign("a", Sign::NonNegative))) << "a has NonNegative";
    EXPECT_FALSE((store.has_sign("a", Sign::Positive))) << "a does not have Positive";
    EXPECT_FALSE((store.has_sign("a", Sign::NonPositive))) << "a does not have NonPositive";
    EXPECT_FALSE((store.has_sign("a", Sign::NonZero))) << "a does not have NonZero";
    EXPECT_FALSE((store.has_sign("a", Sign::Zero))) << "a does not have Zero";
    EXPECT_FALSE((store.has_sign("a", Sign::Negative))) << "a does not have Negative";
}

TEST(PropertyStoreSign, NonpositiveNoExtraImplications) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_sign("b", Sign::NonPositive).has_value())) << "sign declaration succeeds";

    EXPECT_TRUE((store.has_sign("b", Sign::NonPositive))) << "b has NonPositive";
    EXPECT_FALSE((store.has_sign("b", Sign::Negative))) << "b does not have Negative";
    EXPECT_FALSE((store.has_sign("b", Sign::NonNegative))) << "b does not have NonNegative";
    EXPECT_FALSE((store.has_sign("b", Sign::NonZero))) << "b does not have NonZero";
    EXPECT_FALSE((store.has_sign("b", Sign::Zero))) << "b does not have Zero";
    EXPECT_FALSE((store.has_sign("b", Sign::Positive))) << "b does not have Positive";
}

TEST(PropertyStoreSign, NonzeroNoExtraImplications) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_sign("c", Sign::NonZero).has_value())) << "sign declaration succeeds";

    EXPECT_TRUE((store.has_sign("c", Sign::NonZero))) << "c has NonZero";
    EXPECT_FALSE((store.has_sign("c", Sign::Positive))) << "c does not have Positive";
    EXPECT_FALSE((store.has_sign("c", Sign::Negative))) << "c does not have Negative";
    EXPECT_FALSE((store.has_sign("c", Sign::NonNegative))) << "c does not have NonNegative";
    EXPECT_FALSE((store.has_sign("c", Sign::NonPositive))) << "c does not have NonPositive";
    EXPECT_FALSE((store.has_sign("c", Sign::Zero))) << "c does not have Zero";
}

TEST(PropertyStoreSign, IdempotentRedeclaration) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_sign("x", Sign::Positive).has_value())) << "sign declaration succeeds";

    // Re-declaring Positive should be a no-op (no exception)
    EXPECT_TRUE((store.declare_sign("x", Sign::Positive).has_value())) << "sign declaration succeeds";
    EXPECT_TRUE((store.has_sign("x", Sign::Positive))) << "x still has Positive after re-declaration";
    EXPECT_TRUE((store.has_sign("x", Sign::NonNegative))) << "x still has NonNegative after re-declaration";
    EXPECT_TRUE((store.has_sign("x", Sign::NonZero))) << "x still has NonZero after re-declaration";

    // Re-declaring an implied sign should also be a no-op
    EXPECT_TRUE((store.declare_sign("x", Sign::NonNegative).has_value())) << "sign declaration succeeds";
    EXPECT_TRUE((store.has_sign("x", Sign::NonNegative))) << "x still has NonNegative after implied re-declaration";

    EXPECT_TRUE((store.declare_sign("x", Sign::NonZero).has_value())) << "sign declaration succeeds";
    EXPECT_TRUE((store.has_sign("x", Sign::NonZero))) << "x still has NonZero after implied re-declaration";
}

TEST(PropertyStoreSign, ContradictionPositiveNegative) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_sign("x", Sign::Positive).has_value())) << "sign declaration succeeds";

    auto failure_112 = store.declare_sign("x", Sign::Negative);
    EXPECT_TRUE((!failure_112.has_value())) << "Declaring Negative after Positive returns failure";
    // State should be unchanged
    EXPECT_TRUE((store.has_sign("x", Sign::Positive))) << "x still has Positive after failed declaration";
    EXPECT_FALSE((store.has_sign("x", Sign::Negative))) << "x does not have Negative after failed declaration";
}

TEST(PropertyStoreSign, ContradictionPositiveZero) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_sign("x", Sign::Positive).has_value())) << "sign declaration succeeds";

    auto failure_129 = store.declare_sign("x", Sign::Zero);
    EXPECT_TRUE((!failure_129.has_value())) << "Declaring Zero after Positive returns failure";
    EXPECT_FALSE((store.has_sign("x", Sign::Zero))) << "x does not have Zero after failed declaration";
}

TEST(PropertyStoreSign, ContradictionPositiveNonpositive) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_sign("x", Sign::Positive).has_value())) << "sign declaration succeeds";

    auto failure_144 = store.declare_sign("x", Sign::NonPositive);
    EXPECT_TRUE((!failure_144.has_value())) << "Declaring NonPositive after Positive returns failure";
    EXPECT_FALSE((store.has_sign("x", Sign::NonPositive))) << "x does not have NonPositive after failed declaration";
}

TEST(PropertyStoreSign, ContradictionNegativeZero) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_sign("x", Sign::Negative).has_value())) << "sign declaration succeeds";

    auto failure_159 = store.declare_sign("x", Sign::Zero);
    EXPECT_TRUE((!failure_159.has_value())) << "Declaring Zero after Negative returns failure";
    EXPECT_FALSE((store.has_sign("x", Sign::Zero))) << "x does not have Zero after failed declaration";
}

TEST(PropertyStoreSign, ContradictionNegativeNonnegative) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_sign("x", Sign::Negative).has_value())) << "sign declaration succeeds";

    auto failure_174 = store.declare_sign("x", Sign::NonNegative);
    EXPECT_TRUE((!failure_174.has_value())) << "Declaring NonNegative after Negative returns failure";
    EXPECT_FALSE((store.has_sign("x", Sign::NonNegative))) << "x does not have NonNegative after failed declaration";
}

TEST(PropertyStoreSign, ContradictionNonnegativeNegative) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_sign("x", Sign::NonNegative).has_value())) << "sign declaration succeeds";

    auto failure_189 = store.declare_sign("x", Sign::Negative);
    EXPECT_TRUE((!failure_189.has_value())) << "Declaring Negative after NonNegative returns failure";
    EXPECT_FALSE((store.has_sign("x", Sign::Negative))) << "x does not have Negative after failed declaration";
}

TEST(PropertyStoreSign, ContradictionNonpositivePositive) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_sign("x", Sign::NonPositive).has_value())) << "sign declaration succeeds";

    auto failure_204 = store.declare_sign("x", Sign::Positive);
    EXPECT_TRUE((!failure_204.has_value())) << "Declaring Positive after NonPositive returns failure";
    EXPECT_FALSE((store.has_sign("x", Sign::Positive))) << "x does not have Positive after failed declaration";
}

TEST(PropertyStoreSign, ContradictionZeroNonzero) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_sign("x", Sign::Zero).has_value())) << "sign declaration succeeds";

    auto failure_219 = store.declare_sign("x", Sign::NonZero);
    EXPECT_TRUE((!failure_219.has_value())) << "Declaring NonZero after Zero returns failure";
    EXPECT_FALSE((store.has_sign("x", Sign::NonZero))) << "x does not have NonZero after failed declaration";
}

TEST(PropertyStoreSign, ContradictionViaImpliedSigns) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_sign("x", Sign::Positive).has_value())) << "sign declaration succeeds";
    // x now has: Positive, NonNegative, NonZero

    // Declaring Negative should fail because Negative contradicts NonNegative (implied)
    auto failure_237 = store.declare_sign("x", Sign::Negative);
    EXPECT_TRUE((!failure_237.has_value())) << "Declaring Negative contradicts implied NonNegative from Positive";
}

TEST(PropertyStoreSign, ContradictionImpliedAgainstExisting) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_sign("x", Sign::NonZero).has_value())) << "sign declaration succeeds";
    // x now has: NonZero

    // Declaring Zero should fail because Zero contradicts NonZero
    auto failure_253 = store.declare_sign("x", Sign::Zero);
    EXPECT_TRUE((!failure_253.has_value())) << "Declaring Zero contradicts existing NonZero";
}

TEST(PropertyStoreSign, CompatibleSignsCanCoexist) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_sign("x", Sign::NonNegative).has_value())) << "sign declaration succeeds";
    EXPECT_TRUE((store.declare_sign("x", Sign::NonZero).has_value())) << "sign declaration succeeds";

    EXPECT_TRUE((store.has_sign("x", Sign::NonNegative))) << "x has NonNegative";
    EXPECT_TRUE((store.has_sign("x", Sign::NonZero))) << "x has NonZero";
}

TEST(PropertyStoreSign, PositiveAfterNonnegative) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_sign("x", Sign::NonNegative).has_value())) << "sign declaration succeeds";
    EXPECT_TRUE((store.declare_sign("x", Sign::Positive).has_value())) << "sign declaration succeeds";

    EXPECT_TRUE((store.has_sign("x", Sign::Positive))) << "x has Positive";
    EXPECT_TRUE((store.has_sign("x", Sign::NonNegative))) << "x has NonNegative";
    EXPECT_TRUE((store.has_sign("x", Sign::NonZero))) << "x has NonZero (implied by Positive)";
}

TEST(PropertyStoreSign, GetSignsReturnsAllStored) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_sign("x", Sign::Positive).has_value())) << "sign declaration succeeds";

    auto signs = store.get_signs("x");
    EXPECT_TRUE((signs.count(Sign::Positive) > 0)) << "get_signs includes Positive";
    EXPECT_TRUE((signs.count(Sign::NonNegative) > 0)) << "get_signs includes NonNegative";
    EXPECT_TRUE((signs.count(Sign::NonZero) > 0)) << "get_signs includes NonZero";
    EXPECT_TRUE((signs.size() == 3)) << "get_signs has exactly 3 signs for Positive";
}

TEST(PropertyStoreSign, UndeclaredSymbolHasNoSigns) {
    PropertyStore store;

    EXPECT_FALSE((store.has_sign("unknown", Sign::Positive))) << "unknown has no Positive";
    EXPECT_FALSE((store.has_sign("unknown", Sign::Negative))) << "unknown has no Negative";
    EXPECT_FALSE((store.has_sign("unknown", Sign::NonNegative))) << "unknown has no NonNegative";
    EXPECT_FALSE((store.has_sign("unknown", Sign::NonPositive))) << "unknown has no NonPositive";
    EXPECT_FALSE((store.has_sign("unknown", Sign::Zero))) << "unknown has no Zero";
    EXPECT_FALSE((store.has_sign("unknown", Sign::NonZero))) << "unknown has no NonZero";

    auto signs = store.get_signs("unknown");
    EXPECT_TRUE((signs.empty())) << "get_signs returns empty set for undeclared symbol";
}

TEST(PropertyStoreSign, ZeroDomainPromotionDoesNotOverrideMoreSpecific) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_domain("x", Domain::Natural).has_value())) << "domain declaration succeeds";
    EXPECT_TRUE((store.declare_sign("x", Sign::Zero).has_value())) << "sign declaration succeeds";

    // Natural is more specific than Integer, so domain should remain Natural
    EXPECT_TRUE((store.get_domain("x") == Domain::Natural)) << "Domain remains Natural (more specific than Integer)";
    EXPECT_TRUE((store.has_sign("x", Sign::Zero))) << "x has Zero";
}

static void test_checked_property_store_contracts_signs(PropertyStore &store) {
    auto bad_symbol = store.declare_sign_checked("", Sign::Positive);
    EXPECT_TRUE((!bad_symbol.has_value())) << "checked declare_sign rejects empty symbol";
    EXPECT_TRUE((bad_symbol.error().code == CasErrc::InvalidArgument)) << "checked declare_sign reports InvalidArgument for empty symbol";
    EXPECT_FALSE((store.has_sign("", Sign::Positive))) << "failed checked declare_sign does not create empty-symbol fact";

    auto positive = store.declare_sign_checked("x", Sign::Positive);
    EXPECT_TRUE((positive.has_value())) << "checked declare_sign succeeds";
    EXPECT_TRUE((store.has_sign("x", Sign::Positive))) << "checked declare_sign stores declared sign";

    auto contradiction = store.declare_sign_checked("x", Sign::Negative);
    EXPECT_TRUE((!contradiction.has_value())) << "checked declare_sign rejects contradiction";
    EXPECT_TRUE((contradiction.error().code == CasErrc::InvalidArgument)) << "checked declare_sign reports InvalidArgument for contradiction";
    EXPECT_TRUE((store.has_sign("x", Sign::Positive))) << "failed checked declare_sign preserves previous sign";
    EXPECT_FALSE((store.has_sign("x", Sign::Negative))) << "failed checked declare_sign does not apply conflicting sign";
}

static void test_checked_property_store_contracts_domains_and_parity(PropertyStore &store) {
    auto domain = store.declare_domain_checked("n", Domain::Integer);
    EXPECT_TRUE((domain.has_value())) << "checked declare_domain succeeds";
    EXPECT_TRUE((store.has_domain("n", Domain::Integer))) << "checked declare_domain stores domain";

    auto parity = store.declare_parity_checked("n", Parity::Even);
    EXPECT_TRUE((parity.has_value())) << "checked declare_parity succeeds";
    auto parity_conflict = store.declare_parity_checked("n", Parity::Odd);
    EXPECT_TRUE((!parity_conflict.has_value())) << "checked declare_parity rejects contradiction";
    EXPECT_TRUE((parity_conflict.error().code == CasErrc::InvalidArgument)) << "checked declare_parity reports InvalidArgument for contradiction";
    EXPECT_TRUE((store.get_parity("n") == Parity::Even)) << "failed checked declare_parity preserves previous parity";
}

static void test_checked_property_store_contracts_boundedness_and_algebraicity(PropertyStore &store) {
    auto bounded = store.declare_bounded_checked("f", Boundedness::Bounded);
    EXPECT_TRUE((bounded.has_value())) << "checked declare_bounded succeeds";
    auto bounded_conflict = store.declare_bounded_checked("f", Boundedness::Unbounded);
    EXPECT_TRUE((!bounded_conflict.has_value())) << "checked declare_bounded rejects contradiction";
    EXPECT_TRUE((store.get_boundedness("f") == Boundedness::Bounded)) << "failed checked declare_bounded preserves previous boundedness";

    auto transcendental = store.declare_transcendental_checked("t");
    EXPECT_TRUE((transcendental.has_value())) << "checked declare_transcendental succeeds";
    auto algebraic_conflict = store.declare_domain_checked("t", Domain::Algebraic);
    EXPECT_TRUE((!algebraic_conflict.has_value())) << "checked declare_domain rejects transcendental/algebraic conflict";
    EXPECT_TRUE((store.is_transcendental("t"))) << "failed checked declare_domain preserves transcendental marker";
}

static void test_checked_property_store_contracts_finiteness_and_definiteness(PropertyStore &store) {
    auto finite = store.declare_finiteness_checked("seq", Finiteness::Finite);
    EXPECT_TRUE((finite.has_value())) << "checked declare_finiteness succeeds";
    auto divergent = store.declare_finiteness_checked("seq", Finiteness::Divergent);
    EXPECT_TRUE((!divergent.has_value())) << "checked declare_finiteness rejects contradiction";
    EXPECT_TRUE((store.get_finiteness("seq") == Finiteness::Finite)) << "failed checked declare_finiteness preserves previous finiteness";

    auto pd = store.declare_definiteness_checked("A", Definiteness::PositiveDefinite);
    EXPECT_TRUE((pd.has_value())) << "checked declare_definiteness succeeds";
    auto indefinite = store.declare_definiteness_checked("A", Definiteness::Indefinite);
    EXPECT_TRUE((!indefinite.has_value())) << "checked declare_definiteness rejects contradiction";
    EXPECT_TRUE((store.get_definiteness("A") == Definiteness::PositiveDefinite)) << "failed checked declare_definiteness preserves previous definiteness";
}

static void test_checked_property_store_contracts_periodicity(PropertyStore &store) {
    const SymbolicExpr one = *SymbolicExpr::number(1);
    const SymbolicExpr zero = *SymbolicExpr::number(0);
    const auto revision = store.revision();

    auto empty_symbol = store.declare_periodic_checked("", "x", one);
    ASSERT_FALSE(empty_symbol);
    EXPECT_EQ(empty_symbol.error().code, CasErrc::InvalidArgument);
    auto empty_variable = store.declare_periodic_checked("g", "", one);
    ASSERT_FALSE(empty_variable);
    EXPECT_EQ(empty_variable.error().code, CasErrc::InvalidArgument);
    auto nonpositive = store.declare_periodic_checked("g", "x", zero);
    ASSERT_FALSE(nonpositive);
    EXPECT_EQ(nonpositive.error().code, CasErrc::InvalidArgument);
    EXPECT_EQ(store.revision(), revision);
    EXPECT_FALSE(store.is_periodic("g", "x"));

    const SymbolicExpr two = *SymbolicExpr::number(2);
    ASSERT_TRUE(store.declare_periodic_checked("g", "x", two));
    EXPECT_TRUE(store.is_periodic("g", "x"));
    auto stored = store.get_period("g", "x");
    ASSERT_TRUE(stored);
    EXPECT_DOUBLE_EQ(stored->to_numeric(), 2.0);
}

TEST(PropertyStoreSign, CheckedPropertyStoreContracts) {
    PropertyStore store;

    test_checked_property_store_contracts_signs(store);
    test_checked_property_store_contracts_domains_and_parity(store);
    test_checked_property_store_contracts_boundedness_and_algebraicity(store);
    test_checked_property_store_contracts_finiteness_and_definiteness(store);
    test_checked_property_store_contracts_periodicity(store);
}

TEST(PropertyStoreSign, LegacyDeclarationsDelegateTransactionally) {
    PropertyStore store;
    EXPECT_TRUE((store.declare_sign("x", Sign::Positive).has_value())) << "sign declaration succeeds";
    const auto symbols_before = store.get_all_symbols();

    auto contradiction = store.declare_sign("x", Sign::Negative);
    EXPECT_TRUE((!contradiction.has_value())) << "canonical sign declaration returns contradiction";
    EXPECT_TRUE((contradiction.error().code == CasErrc::InvalidArgument)) << "sign contradiction reports InvalidArgument";
    EXPECT_TRUE((store.has_sign("x", Sign::Positive))) << "canonical sign failure preserves the established sign";
    EXPECT_FALSE((store.has_sign("x", Sign::Negative))) << "canonical sign failure does not commit a contradictory sign";

    auto empty_symbol = store.declare_domain("", Domain::Real);
    EXPECT_TRUE((!empty_symbol.has_value())) << "canonical domain declaration rejects an empty symbol";
    EXPECT_TRUE((empty_symbol.error().code == CasErrc::InvalidArgument)) << "empty domain symbol reports InvalidArgument";
    EXPECT_TRUE((store.get_all_symbols() == symbols_before)) << "failed canonical declaration does not create an empty property record";
}
