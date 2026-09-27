#pragma once

#include "property_store.hpp"

namespace LMCAS::property_detail {

std::string domain_str(Domain domain);
PropertyStoreResult invalid_empty_symbol(const std::string& operation);

template <typename F>
PropertyStoreResult checked_property_update(PropertyStore& store,
                                            const std::string& symbol,
                                            const std::string& operation,
                                            F&& update) {
    if (symbol.empty()) {
        return invalid_empty_symbol(operation);
    }
    try {
        PropertyStore candidate = store;
        update(candidate);
        store = std::move(candidate);
    } catch (const std::bad_alloc&) {
        return PropertyStoreResult::failure(
            CasErrc::ResourceLimit, "property-store allocation failed", operation);
    } catch (const std::invalid_argument& ex) {
        return PropertyStoreResult::failure(CasErrc::InvalidArgument, ex.what(), operation);
    } catch (const std::exception& ex) {
        return PropertyStoreResult::failure(CasErrc::InternalInvariant, ex.what(), operation);
    }
    return PropertyStoreResult::success();
}

}
