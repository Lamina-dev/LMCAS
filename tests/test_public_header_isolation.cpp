#ifndef LMCAS_HEADER_UNDER_TEST
#error "LMCAS_HEADER_UNDER_TEST must name one public header"
#endif

#include LMCAS_HEADER_UNDER_TEST

#ifdef LMCAS_INTERNAL_AST_INCLUDED
#error "public module headers must not include the AST implementation"
#endif

#include <type_traits>

using namespace LMCAS;

template <typename T, typename = void>
struct HeaderTypeIsComplete : std::false_type {};

template <typename T>
struct HeaderTypeIsComplete<T, std::void_t<decltype(sizeof(T))>> : std::true_type {};

static_assert(HeaderTypeIsComplete<SymbolicExpr>::value,
              "public module headers must expose the stable expression type");

int lmcas_public_header_isolation_anchor() {
    return 0;
}
