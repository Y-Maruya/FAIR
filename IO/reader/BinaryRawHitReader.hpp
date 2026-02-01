#pragma once 
#include "common/edm/RawHit.hpp"
#include "common/edm/RawData.hpp"

class BinaryRawHitReader{
public:
    BinaryRawHitReader(std::string filename);
    ~BinaryRawHitReader();

    bool next(std::vector<AHCALRawHit>& out_hits, AHCALTLURawData& out_tlu);

    long long entry() const { return m_entry; }

private:
    std::unique_ptr<std::ifstream> m_ifs;
    long long m_entry;
    bool m_eof_good = true;
    bool has_tlu = false;
    bool has_ahcal = false;
    std::vector<int> last_timestamp = std::vector<int>(6, -1);
};