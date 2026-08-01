#pragma once
// The seam between a bitstream parser and VK_KHR_video_decode_*.
//
// A hardware video decoder does the entropy decoding, motion compensation and
// filtering; what it does NOT do is track the reference-picture state, and that
// is what a `CodecDecoder` is for. Per coded frame it produces
//
//   * the bitstream the driver should see (Annex-B for H.264/H.265, OBUs for
//     AV1), plus per-slice offsets into it,
//   * the codec-specific picture-info structure that hangs off
//     VkVideoDecodeInfoKHR,
//   * which DPB slot the picture activates and which slots it references.
//
// DPB slot allocation is codec business (H.264 sliding-window/MMCO marking,
// H.265 reference picture sets, AV1 refresh_frame_flags), so it lives behind
// this interface and VulkanVideo.cpp only ever maps a slot index to an image.

#include "video/Demuxer.h"

#include <vulkan/vulkan.h>

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace video {

struct StreamFormat {
    int width = 0, height = 0;              // display size, after cropping
    int coded_width = 0, coded_height = 0;  // what the session must accommodate
    int bit_depth = 8;
    int chroma_format = 1;                  // 0 mono, 1 4:2:0, 2 4:2:2, 3 4:4:4
    uint32_t max_dpb_slots = 0;             // including the picture being decoded
    uint32_t max_active_references = 0;
    uint32_t max_reorder = 0;               // pictures the display queue may hold
    bool full_range = false;                // VUI video_full_range_flag
    int  matrix_coefficients = 2;           // H.273; 2 == unspecified
    const char* profile_name = "";
};

// One decoded picture's worth of driver-facing state. All pointers refer to
// storage owned by the CodecDecoder and stay valid until the next decodeFrame().
struct PictureInfo {
    const void* decode_pnext = nullptr;   // VkVideoDecodeH264PictureInfoKHR, ...

    struct Ref {
        int32_t     slot = -1;
        const void* std_ref = nullptr;    // StdVideoDecodeH264ReferenceInfo, ...
    };
    int32_t     setup_slot = -1;          // DPB slot this picture activates, or -1
    const void* setup_std_ref = nullptr;
    std::vector<Ref> refs;                // active references, for pReferenceSlots

    int64_t poc = 0;                      // display order within the sequence
    bool    output = true;                // false for an AV1 no-show frame
    int32_t show_existing_slot = -1;      // AV1 show_existing_frame -> re-output
    bool    params_changed = false;       // in-band parameter sets updated
    bool    sequence_restart = false;     // IDR/IRAP: flush the reorder queue
    // An AV1 temporal unit routinely carries several coded frames -- one shown
    // plus the hidden alternate references it was predicted from. When this is
    // set, call decodeFrame() again with a null pointer to get the next one.
    bool    more_in_packet = false;
    // Which DPB slots hold live references AFTER this picture is decoded. The
    // decoder pins the images backing them.
    std::vector<int32_t> live_slots;
};

class CodecDecoder {
public:
    virtual ~CodecDecoder() = default;

    virtual VkVideoCodecOperationFlagBitsKHR operation() const = 0;
    // pNext for VkVideoProfileInfoKHR (VkVideoDecodeH264ProfileInfoKHR, ...).
    virtual const void* profileExt() const = 0;
    virtual const StreamFormat& format() const = 0;

    // Parses the container's codec configuration record. Must be enough to
    // create the video session; parameter sets found in-band later are added
    // through updateParameters().
    virtual bool init(const TrackInfo& track, std::string& error) = 0;

    virtual VkVideoSessionParametersKHR createParameters(VkDevice device,
                                                         VkVideoSessionKHR session,
                                                         std::string& error) = 0;

    // Parses one coded frame. `bitstream` is cleared and refilled with the
    // bytes to hand the driver; `slice_offsets` (H.264/H.265) indexes the start
    // of each slice NAL unit within it. `data == nullptr` continues within the
    // packet last passed in (see PictureInfo::more_in_packet).
    virtual bool decodeFrame(const uint8_t* data, size_t size, int nal_length_size,
                             std::vector<uint8_t>& bitstream,
                             std::vector<uint32_t>& slice_offsets, PictureInfo& out,
                             std::string& error) = 0;

    // Commits the reference marking decided by decodeFrame(). Split out so a
    // frame that fails to record leaves the DPB untouched.
    virtual void commitFrame() = 0;

    virtual void flush() = 0;
};

std::unique_ptr<CodecDecoder> make_codec_decoder(Codec codec);

// Splits length-prefixed NAL units into Annex-B, appending 4-byte start codes.
// Shared by the H.264 and H.265 paths.
void annexb_append(std::vector<uint8_t>& out, const uint8_t* nal, size_t size);

}  // namespace video
