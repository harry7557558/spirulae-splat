#pragma once
// Matroska / WebM (mkv, webm).
//
// Streaming rather than table-driven: the header elements (Info, Tracks) are
// read up front, then clusters are walked lazily. A 500 MB WebM is never held
// in memory, and no index is required -- Matroska stores frames in decode order
// with a presentation timestamp on each block, which is all the reorder queue
// needs.

#include "video/Demuxer.h"

#include <fstream>

namespace video {

class MkvDemuxer : public Demuxer {
public:
    bool open(const std::string& path, std::string& error);

    const std::vector<TrackInfo>& tracks() const override { return infos_; }
    bool selectTrack(int index, std::string& error) override;
    bool next(Packet& out, std::string& error) override;

private:
    struct Track {
        uint64_t  number = 0;
        TrackInfo info;
    };

    // --- primitive EBML reads, all at the current file offset ---
    bool     readId(uint32_t& id);
    bool     readSize(uint64_t& size, bool& unknown);
    uint64_t readUInt(uint64_t len);
    double   readFloat(uint64_t len);
    bool     readBytes(void* dst, uint64_t len);
    void     skip(uint64_t len) { seek(pos_ + len); }
    void     seek(uint64_t off);

    bool parseTracks(uint64_t end, std::string& error);
    bool parseTrackEntry(uint64_t end);
    bool parseInfo(uint64_t end);
    // Advances to the next cluster with payload, setting cluster_ts_/cluster_end_.
    bool nextCluster(std::string& error);

    std::string   path_;
    std::ifstream file_;
    uint64_t file_size_ = 0;
    uint64_t pos_ = 0;

    uint64_t segment_end_ = 0;
    uint64_t first_cluster_ = 0;
    double   timestamp_scale_ = 1e-6;   // seconds per timestamp unit
    double   duration_ = 0.0;           // seconds

    std::vector<Track>     tracks_;
    std::vector<TrackInfo> infos_;
    int      selected_ = -1;

    bool     in_cluster_ = false;
    uint64_t cluster_end_ = 0;
    int64_t  cluster_ts_ = 0;
    int64_t  packet_index_ = 0;
};

}  // namespace video
