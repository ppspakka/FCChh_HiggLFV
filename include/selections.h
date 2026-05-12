#pragma once

#include "iselection.h"

namespace hlfv {

// Concrete selections

struct EmptySelection : public ISelection {
    std::string name() const override;
    bool apply(const Event& evt, Meta& meta, const Parameters& cfg) override;
};

struct FinalState_NoCut : public ISelection {
    std::string name() const override;
    bool apply(const Event& evt, Meta& meta, const Parameters& cfg) override;
};
struct LeptonSelection : public ISelection {
    std::string name() const override;
    bool apply(const Event& evt, Meta& meta, const Parameters& cfg) override;
};
struct JetVeto : public ISelection {
    std::string name() const override;
    bool apply(const Event& evt, Meta& meta, const Parameters& cfg) override;
};
struct JetCategorization : public ISelection {
    std::string name() const override;
    bool apply(const Event& evt, Meta& meta, const Parameters& cfg) override;
};
struct LeptonPTCut : public ISelection {
    std::string name() const override;
    bool apply(const Event& evt, Meta& meta, const Parameters& cfg) override;
};
struct DeltaPhiMETLepton : public ISelection {
    std::string name() const override;
    bool apply(const Event& evt, Meta& meta, const Parameters& cfg) override;
};


} // namespace hlfv
