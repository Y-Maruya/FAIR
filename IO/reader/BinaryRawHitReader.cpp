#include "BinaryRawHitReader.hpp"

#include "EventFormats/AHCALDataFragment.hpp"
#include "EventFormats/DAQFormats.hpp"
#include "common/Logger.hpp"
#include <iostream>
#include <cstdint>
typedef unsigned char uchar_t;
namespace tlu{
  class fmctludata{
  public:
    fmctludata(uint64_t wl, uint64_t wh, uint64_t we):  // wl -> wh
      eventtype((wl>>60)&0xf),
      input0((wl>>48)&0x1),
      input1((wl>>49)&0x1),
      input2((wl>>50)&0x1),
      input3((wl>>51)&0x1),
      input4((wl>>52)&0x1),
      input5((wl>>53)&0x1),
      timestamp(wl&0xffffffffffff),
      sc0((wh>>56)&0xff),
      sc1((wh>>48)&0xff),
      sc2((wh>>40)&0xff),
      sc3((wh>>32)&0xff),
      sc4((we>>56)&0xff),
      sc5((we>>48)&0xff),
      eventnumber(wh&0xffffffff){
    }

    fmctludata(uint32_t w0, uint32_t w1, uint32_t w2, uint32_t w3, uint32_t w4, uint32_t w5): // w0 w1 w2 w3  wl= w0 w1; wh= w2 w3
       eventtype((w0>>28)&0xf),
       input0((w0>>16)&0x1),
       input1((w0>>17)&0x1),
       input2((w0>>18)&0x1),
       input3((w0>>19)&0x1),
       input4((w0>>20)&0x1),
       input5((w0>>21)&0x1),
       timestamp(((uint64_t(w0&0x0000ffff))<<32) + w1),
       // timestamp(w1),
       sc0((w2>>24)&0xff),
       sc1((w2>>16)&0xff),
       sc2((w2>>8)&0xff),
       sc3(w2&0xff),
       sc4((w4>>24)&0xff),
       sc5((w4>>16)&0xff),
       eventnumber(w3),
       timestamp1(w0&0xffff){
        (void)w5; // suppress unused-parameter
        payload[0]= w0;
        payload[1]= w1;
        payload[2]= w2;
        payload[3]= w3;
        payload[4]= w4;
    }

    uint32_t payload[5];

    uchar_t eventtype;
    uchar_t input0;
    uchar_t input1;
    uchar_t input2;
    uchar_t input3;
    uchar_t input4;
    uchar_t input5;
    uint64_t timestamp;
    uchar_t sc0;
    uchar_t sc1;
    uchar_t sc2;
    uchar_t sc3;
    uchar_t sc4;
    uchar_t sc5;
    uint32_t eventnumber; // trigger number
    uint64_t timestamp1;

  };

  std::ostream &operator<<(std::ostream &s, fmctludata &d);
}

AHCALRawHit convert(int index, const AHCALHit& hit){
    AHCALRawHit rawhit;
    rawhit.cellID = hit.layer*100000 + hit.asic*10000 + hit.channel;
    rawhit.hg_adc = hit.hg_adc;
    rawhit.lg_adc = hit.lg_adc;
    rawhit.hittag = hit.hittag;
    rawhit.bcid = hit.bcid;
    rawhit.index = index;
    return rawhit;
}

BinaryRawHitReader::BinaryRawHitReader(std::string filename)
{
    m_ifs = std::make_unique<std::ifstream>(filename,  std::ios::binary);
    if(!m_ifs->is_open()) {
        LOG_ERROR("Failed to open file: {}", filename);
        throw std::runtime_error("Failed to open file");
    }
    LOG_INFO("Opened binary raw hit file: {}", filename);
}

BinaryRawHitReader::~BinaryRawHitReader()
{
    if(m_ifs && m_ifs->is_open()) {
        m_ifs->close();
        LOG_INFO("Closed binary raw hit file.");
    }
}

bool BinaryRawHitReader::next(std::vector<AHCALRawHit>& out_hits, AHCALTLURawData& out_tlu) {
    ++m_entry;
    m_eof_good = m_ifs->good() and m_ifs->peek()!=EOF;
    if (!m_eof_good) return false;
    out_hits.clear();
    out_tlu = AHCALTLURawData();
    try {
        DAQFormats::EventFull event(*m_ifs);
        out_tlu.RunNo = event.run_number();
        // TLUDATA
        for (const auto id: event.getFragmentIDs()){
            const DAQFormats::EventFragment* frag = event.find_fragment(id);
            if (frag->source_id() == DAQFormats::SourceIDs::TLUSourceID){
                const uint8_t *payload_8bit = frag->payload<const uint8_t*>();
                if (payload_8bit == nullptr) {
                    LOG_WARN("TLU Fragment payload is null.");
                    return false;
                }
                const uint32_t* payload = reinterpret_cast<const uint32_t*>(payload_8bit);
                if (frag->payload_size() != 20) {
                    LOG_WARN("Unexpected payload size: {} bytes. Expected 20 bytes.", frag->payload_size());
                    return false;
                }
                
                // Stack allocation instead of new/delete
                tlu::fmctludata data(payload[0], payload[1], payload[2], payload[3], payload[4], 0);
                out_tlu.TriggerID = data.eventnumber;
                out_tlu.Timestamp = data.timestamp;
                out_tlu.BCID_TLU = static_cast<int>(out_tlu.Timestamp & 0xFFFF);
                
                out_tlu.Inputs.clear();
                out_tlu.Inputs.reserve(6);
                out_tlu.Inputs.push_back(static_cast<int>(data.input0));
                out_tlu.Inputs.push_back(static_cast<int>(data.input1));
                out_tlu.Inputs.push_back(static_cast<int>(data.input2));
                out_tlu.Inputs.push_back(static_cast<int>(data.input3));
                out_tlu.Inputs.push_back(static_cast<int>(data.input4));
                out_tlu.Inputs.push_back(static_cast<int>(data.input5));
                
                out_tlu.FineTimestamps.clear();
                out_tlu.FineTimestamps.reserve(6);
                std::vector<int> current_timestamps(6);
                current_timestamps[0] = static_cast<int>(data.sc0);
                current_timestamps[1] = static_cast<int>(data.sc1);
                current_timestamps[2] = static_cast<int>(data.sc2);
                current_timestamps[3] = static_cast<int>(data.sc3);
                current_timestamps[4] = static_cast<int>(data.sc4);
                current_timestamps[5] = static_cast<int>(data.sc5);
                
                for (int i = 0; i < 6; ++i) {
                    if (current_timestamps[i] == last_timestamp[i] && out_tlu.Inputs[i] == 0) {
                        out_tlu.FineTimestamps.push_back(-1);
                    } else {
                        out_tlu.FineTimestamps.push_back(current_timestamps[i]);
                        last_timestamp[i] = current_timestamps[i];
                    }
                }
                has_tlu = true;
            }else if (frag->source_id() == DAQFormats::SourceIDs::AHCALSourceID){
                const uint8_t* payload_8bit = frag->payload<const uint8_t*>();
                if (payload_8bit == nullptr) {
                    LOG_WARN("AHCAL Fragment payload is null.");
                    return false;
                }
                const std::vector<uint8_t> payload(payload_8bit, payload_8bit + frag->payload_size());
                AHCALDataFragment ahcalFrag(payload);
                const auto& hits = ahcalFrag.hits();
                out_tlu.CycleID = static_cast<int>(ahcalFrag.cycleID());
                out_tlu.Event_Time = static_cast<int>(ahcalFrag.timestamp());
                int i = 0;
                for (AHCALHit hit : hits){
                    out_hits.push_back(convert(i, hit));
                    ++i;
                }
                has_ahcal = true; 
            }                
        }
    } catch (const std::exception& e) {
        LOG_ERROR("Error reading event at entry {}: {}", m_entry, e.what());
        return false;
    }
    return has_tlu || has_ahcal;
}
