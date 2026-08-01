// H.264 / AVC picture-parameter driver for VK_KHR_video_decode_h264.
//
// The driver entropy-decodes the slices, so this file never touches macroblock
// data. What it must get exactly right is everything the driver cannot recover
// on its own: the sequence and picture parameter sets, the picture order count,
// and the reference-picture marking that decides which DPB slot holds what.
//
// Frame pictures only. Field and MBAFF-with-field-pictures streams are rejected
// by name rather than decoded into interleaved garbage -- Vulkan can express
// them, but nothing that feeds this library produces them.

#include "nn/core/Log.h"
#include "video/BitReader.h"
#include "video/CodecDecoder.h"
#include "video/VideoApi.h"

#include <algorithm>
#include <array>
#include <cstring>

namespace video {

namespace {

constexpr int kMaxDpb = 17;   // 16 references + the current picture

// Table A-1: MaxDpbMbs, indexed by StdVideoH264LevelIdc.
uint32_t max_dpb_mbs(StdVideoH264LevelIdc level) {
    switch (level) {
        case STD_VIDEO_H264_LEVEL_IDC_1_0: return 396;
        case STD_VIDEO_H264_LEVEL_IDC_1_1: return 900;
        case STD_VIDEO_H264_LEVEL_IDC_1_2:
        case STD_VIDEO_H264_LEVEL_IDC_1_3:
        case STD_VIDEO_H264_LEVEL_IDC_2_0: return 2376;
        case STD_VIDEO_H264_LEVEL_IDC_2_1: return 4752;
        case STD_VIDEO_H264_LEVEL_IDC_2_2:
        case STD_VIDEO_H264_LEVEL_IDC_3_0: return 8100;
        case STD_VIDEO_H264_LEVEL_IDC_3_1: return 18000;
        case STD_VIDEO_H264_LEVEL_IDC_3_2: return 20480;
        case STD_VIDEO_H264_LEVEL_IDC_4_0:
        case STD_VIDEO_H264_LEVEL_IDC_4_1: return 32768;
        case STD_VIDEO_H264_LEVEL_IDC_4_2: return 34816;
        case STD_VIDEO_H264_LEVEL_IDC_5_0: return 110400;
        case STD_VIDEO_H264_LEVEL_IDC_5_1:
        case STD_VIDEO_H264_LEVEL_IDC_5_2: return 184320;
        default: return 696320;
    }
}

StdVideoH264LevelIdc level_from_idc(uint32_t level_idc, bool constraint_set3) {
    switch (level_idc) {
        case 10: return STD_VIDEO_H264_LEVEL_IDC_1_0;
        case 11: return constraint_set3 ? STD_VIDEO_H264_LEVEL_IDC_1_0
                                        : STD_VIDEO_H264_LEVEL_IDC_1_1;
        case 12: return STD_VIDEO_H264_LEVEL_IDC_1_2;
        case 13: return STD_VIDEO_H264_LEVEL_IDC_1_3;
        case 20: return STD_VIDEO_H264_LEVEL_IDC_2_0;
        case 21: return STD_VIDEO_H264_LEVEL_IDC_2_1;
        case 22: return STD_VIDEO_H264_LEVEL_IDC_2_2;
        case 30: return STD_VIDEO_H264_LEVEL_IDC_3_0;
        case 31: return STD_VIDEO_H264_LEVEL_IDC_3_1;
        case 32: return STD_VIDEO_H264_LEVEL_IDC_3_2;
        case 40: return STD_VIDEO_H264_LEVEL_IDC_4_0;
        case 41: return STD_VIDEO_H264_LEVEL_IDC_4_1;
        case 42: return STD_VIDEO_H264_LEVEL_IDC_4_2;
        case 50: return STD_VIDEO_H264_LEVEL_IDC_5_0;
        case 51: return STD_VIDEO_H264_LEVEL_IDC_5_1;
        case 52: return STD_VIDEO_H264_LEVEL_IDC_5_2;
        case 60: return STD_VIDEO_H264_LEVEL_IDC_6_0;
        case 61: return STD_VIDEO_H264_LEVEL_IDC_6_1;
        default: return STD_VIDEO_H264_LEVEL_IDC_6_2;
    }
}

// 7.3.2.1.1.1 scaling_list()
void parse_scaling_list(BitReader& br, uint8_t* list, int size, bool& use_default) {
    int last = 8, next = 8;
    use_default = false;
    for (int j = 0; j < size; ++j) {
        if (next != 0) {
            const int delta = br.se();
            next = (last + delta + 256) % 256;
            if (j == 0 && next == 0) use_default = true;
        }
        list[j] = (uint8_t)(next == 0 ? last : next);
        last = list[j];
    }
}

struct Sps {
    StdVideoH264SequenceParameterSet    std{};
    StdVideoH264ScalingLists            scaling{};
    StdVideoH264SequenceParameterSetVui vui{};
    std::vector<int32_t>                offset_for_ref_frame;
    // Derived once, used everywhere.
    uint32_t max_frame_num = 16;
    uint32_t max_poc_lsb = 16;
    int      width = 0, height = 0;
    uint32_t max_dpb_frames = 16;
    uint32_t max_num_reorder = 16;
    bool     valid = false;
};

struct Pps {
    StdVideoH264PictureParameterSet std{};
    StdVideoH264ScalingLists        scaling{};
    bool valid = false;
};

struct SliceHeader {
    uint32_t first_mb = 0;
    uint32_t slice_type = 0;      // 0..4 after the %5
    uint32_t pps_id = 0;
    uint32_t frame_num = 0;
    bool     field_pic = false;
    bool     bottom_field = false;
    uint32_t idr_pic_id = 0;
    uint32_t poc_lsb = 0;
    int32_t  delta_poc_bottom = 0;
    int32_t  delta_poc[2] = {0, 0};
    // dec_ref_pic_marking
    bool     no_output_of_prior_pics = false;
    bool     long_term_reference_flag = false;
    bool     adaptive_marking = false;
    struct Mmco { uint32_t op; uint32_t arg0; uint32_t arg1; };
    std::vector<Mmco> mmco;
};

// A frame held in the decoded picture buffer.
struct DpbFrame {
    bool used = false;
    bool long_term = false;
    int  frame_num = 0;
    int  frame_num_wrap = 0;
    int  pic_num = 0;
    int  long_term_frame_idx = 0;
    int  poc = 0;
    StdVideoDecodeH264ReferenceInfo std_ref{};
};

class H264Decoder final : public CodecDecoder {
public:
    VkVideoCodecOperationFlagBitsKHR operation() const override {
        return VK_VIDEO_CODEC_OPERATION_DECODE_H264_BIT_KHR;
    }
    const void* profileExt() const override { return &profile_; }
    const StreamFormat& format() const override { return format_; }

    bool init(const TrackInfo& track, std::string& error) override;
    VkVideoSessionParametersKHR createParameters(VkDevice device, VkVideoSessionKHR session,
                                                  std::string& error) override;
    bool decodeFrame(const uint8_t* data, size_t size, int nal_length_size,
                     std::vector<uint8_t>& bitstream, std::vector<uint32_t>& slice_offsets,
                     PictureInfo& out, std::string& error) override;
    void commitFrame() override;
    void flush() override;

private:
    bool parseSps(const uint8_t* nal, size_t size, std::string& error);
    bool parsePps(const uint8_t* nal, size_t size, std::string& error);
    bool parseSliceHeader(const uint8_t* nal, size_t size, int nal_ref_idc, bool idr,
                          SliceHeader& sh, std::string& error);
    void updateFormat(const Sps& sps);
    void computePoc(const SliceHeader& sh, const Sps& sps, int nal_ref_idc, bool idr);
    void applyMarking(const SliceHeader& sh, const Sps& sps);
    int  allocSlot();

    std::array<std::unique_ptr<Sps>, 32>  sps_;
    std::array<std::unique_ptr<Pps>, 256> pps_;
    std::array<DpbFrame, kMaxDpb>         dpb_{};

    VkVideoDecodeH264ProfileInfoKHR profile_{
        VK_STRUCTURE_TYPE_VIDEO_DECODE_H264_PROFILE_INFO_KHR, nullptr,
        STD_VIDEO_H264_PROFILE_IDC_HIGH,
        VK_VIDEO_DECODE_H264_PICTURE_LAYOUT_PROGRESSIVE_KHR};
    StreamFormat format_{};

    // Per-picture output, kept alive for the caller's record.
    StdVideoDecodeH264PictureInfo   std_pic_{};
    VkVideoDecodeH264PictureInfoKHR vk_pic_{};
    StdVideoDecodeH264ReferenceInfo setup_ref_{};

    // POC state.
    int  prev_poc_msb_ = 0, prev_poc_lsb_ = 0;
    int  prev_frame_num_ = 0, prev_frame_num_offset_ = 0;
    bool prev_mmco5_ = false;
    int  top_poc_ = 0, bottom_poc_ = 0;

    // Pending marking for the picture just parsed.
    SliceHeader pending_sh_{};
    const Sps*  pending_sps_ = nullptr;
    int         pending_slot_ = -1;
    bool        pending_is_ref_ = false;
    bool        pending_idr_ = false;
    bool        params_dirty_ = false;
    bool        have_sps_ = false;
};

// ---------------------------------------------------------------------------
// Parameter sets
// ---------------------------------------------------------------------------

bool H264Decoder::parseSps(const uint8_t* nal, size_t size, std::string& error) {
    RbspReader br(nal, size, 1);
    auto sps = std::make_unique<Sps>();
    StdVideoH264SequenceParameterSet& s = sps->std;

    const uint32_t profile_idc = br.u(8);
    const uint32_t constraint = br.u(8);
    const uint32_t level_idc = br.u(8);
    s.profile_idc = (StdVideoH264ProfileIdc)profile_idc;
    s.flags.constraint_set0_flag = (constraint >> 7) & 1;
    s.flags.constraint_set1_flag = (constraint >> 6) & 1;
    s.flags.constraint_set2_flag = (constraint >> 5) & 1;
    s.flags.constraint_set3_flag = (constraint >> 4) & 1;
    s.flags.constraint_set4_flag = (constraint >> 3) & 1;
    s.flags.constraint_set5_flag = (constraint >> 2) & 1;
    s.level_idc = level_from_idc(level_idc, s.flags.constraint_set3_flag != 0);

    const uint32_t sps_id = br.ue();
    if (sps_id >= 32) {
        error = "H.264 SPS id " + std::to_string(sps_id) + " out of range";
        return false;
    }
    s.seq_parameter_set_id = (uint8_t)sps_id;
    s.chroma_format_idc = STD_VIDEO_H264_CHROMA_FORMAT_IDC_420;

    if (profile_idc == 100 || profile_idc == 110 || profile_idc == 122 ||
        profile_idc == 244 || profile_idc == 44 || profile_idc == 83 ||
        profile_idc == 86 || profile_idc == 118 || profile_idc == 128 ||
        profile_idc == 138 || profile_idc == 139 || profile_idc == 134 ||
        profile_idc == 135) {
        const uint32_t chroma = br.ue();
        s.chroma_format_idc = (StdVideoH264ChromaFormatIdc)chroma;
        if (chroma == 3) s.flags.separate_colour_plane_flag = br.flag();
        s.bit_depth_luma_minus8 = (uint8_t)br.ue();
        s.bit_depth_chroma_minus8 = (uint8_t)br.ue();
        s.flags.qpprime_y_zero_transform_bypass_flag = br.flag();
        s.flags.seq_scaling_matrix_present_flag = br.flag();
        if (s.flags.seq_scaling_matrix_present_flag) {
            const int n_lists = (chroma != 3) ? 8 : 12;
            for (int i = 0; i < n_lists; ++i) {
                if (!br.flag()) continue;
                bool use_default = false;
                if (i < 6) {
                    parse_scaling_list(br, sps->scaling.ScalingList4x4[i], 16, use_default);
                    sps->scaling.scaling_list_present_mask |= (uint16_t)(1u << i);
                    if (use_default)
                        sps->scaling.use_default_scaling_matrix_mask |= (uint16_t)(1u << i);
                } else {
                    parse_scaling_list(br, sps->scaling.ScalingList8x8[i - 6], 64, use_default);
                    sps->scaling.scaling_list_present_mask |= (uint16_t)(1u << i);
                    if (use_default)
                        sps->scaling.use_default_scaling_matrix_mask |= (uint16_t)(1u << i);
                }
            }
            s.pScalingLists = &sps->scaling;
        }
    }

    s.log2_max_frame_num_minus4 = (uint8_t)br.ue();
    const uint32_t poc_type = br.ue();
    s.pic_order_cnt_type = (StdVideoH264PocType)poc_type;
    if (poc_type == 0) {
        s.log2_max_pic_order_cnt_lsb_minus4 = (uint8_t)br.ue();
    } else if (poc_type == 1) {
        s.flags.delta_pic_order_always_zero_flag = br.flag();
        s.offset_for_non_ref_pic = br.se();
        s.offset_for_top_to_bottom_field = br.se();
        const uint32_t n = br.ue();
        s.num_ref_frames_in_pic_order_cnt_cycle = (uint8_t)n;
        sps->offset_for_ref_frame.resize(n);
        for (uint32_t i = 0; i < n; ++i) sps->offset_for_ref_frame[i] = br.se();
        s.pOffsetForRefFrame = sps->offset_for_ref_frame.data();
    }
    s.max_num_ref_frames = (uint8_t)br.ue();
    s.flags.gaps_in_frame_num_value_allowed_flag = br.flag();
    s.pic_width_in_mbs_minus1 = br.ue();
    s.pic_height_in_map_units_minus1 = br.ue();
    s.flags.frame_mbs_only_flag = br.flag();
    if (!s.flags.frame_mbs_only_flag) s.flags.mb_adaptive_frame_field_flag = br.flag();
    s.flags.direct_8x8_inference_flag = br.flag();
    s.flags.frame_cropping_flag = br.flag();
    if (s.flags.frame_cropping_flag) {
        s.frame_crop_left_offset = br.ue();
        s.frame_crop_right_offset = br.ue();
        s.frame_crop_top_offset = br.ue();
        s.frame_crop_bottom_offset = br.ue();
    }
    s.flags.vui_parameters_present_flag = br.flag();

    sps->max_num_reorder = 16;
    if (s.flags.vui_parameters_present_flag) {
        StdVideoH264SequenceParameterSetVui& v = sps->vui;
        v.flags.aspect_ratio_info_present_flag = br.flag();
        if (v.flags.aspect_ratio_info_present_flag) {
            v.aspect_ratio_idc = (StdVideoH264AspectRatioIdc)br.u(8);
            if (v.aspect_ratio_idc == STD_VIDEO_H264_ASPECT_RATIO_IDC_EXTENDED_SAR) {
                v.sar_width = (uint16_t)br.u(16);
                v.sar_height = (uint16_t)br.u(16);
            }
        }
        v.flags.overscan_info_present_flag = br.flag();
        if (v.flags.overscan_info_present_flag) v.flags.overscan_appropriate_flag = br.flag();
        v.flags.video_signal_type_present_flag = br.flag();
        if (v.flags.video_signal_type_present_flag) {
            v.video_format = (uint8_t)br.u(3);
            v.flags.video_full_range_flag = br.flag();
            v.flags.color_description_present_flag = br.flag();
            if (v.flags.color_description_present_flag) {
                v.colour_primaries = (uint8_t)br.u(8);
                v.transfer_characteristics = (uint8_t)br.u(8);
                v.matrix_coefficients = (uint8_t)br.u(8);
            }
        }
        v.flags.chroma_loc_info_present_flag = br.flag();
        if (v.flags.chroma_loc_info_present_flag) {
            v.chroma_sample_loc_type_top_field = (uint8_t)br.ue();
            v.chroma_sample_loc_type_bottom_field = (uint8_t)br.ue();
        }
        v.flags.timing_info_present_flag = br.flag();
        if (v.flags.timing_info_present_flag) {
            v.num_units_in_tick = br.u(32);
            v.time_scale = br.u(32);
            v.flags.fixed_frame_rate_flag = br.flag();
        }
        // HRD parameters are parsed only far enough to stay in sync; the
        // decoder has no use for them.
        auto skip_hrd = [&br]() {
            const uint32_t cpb_cnt = br.ue() + 1;
            br.u(4);
            br.u(4);
            for (uint32_t i = 0; i < cpb_cnt && i < 32; ++i) {
                br.ue();
                br.ue();
                br.bit();
            }
            br.u(5);
            br.u(5);
            br.u(5);
            br.u(5);
        };
        v.flags.nal_hrd_parameters_present_flag = br.flag();
        if (v.flags.nal_hrd_parameters_present_flag) skip_hrd();
        v.flags.vcl_hrd_parameters_present_flag = br.flag();
        if (v.flags.vcl_hrd_parameters_present_flag) skip_hrd();
        if (v.flags.nal_hrd_parameters_present_flag || v.flags.vcl_hrd_parameters_present_flag)
            br.bit();  // low_delay_hrd_flag
        br.bit();      // pic_struct_present_flag
        v.flags.bitstream_restriction_flag = br.flag();
        if (v.flags.bitstream_restriction_flag) {
            br.bit();  // motion_vectors_over_pic_boundaries_flag
            br.ue();   // max_bytes_per_pic_denom
            br.ue();   // max_bits_per_mb_denom
            br.ue();   // log2_max_mv_length_horizontal
            br.ue();   // log2_max_mv_length_vertical
            v.max_num_reorder_frames = (uint8_t)br.ue();
            v.max_dec_frame_buffering = (uint8_t)br.ue();
            sps->max_num_reorder = v.max_num_reorder_frames;
        }
        s.pSequenceParameterSetVui = &sps->vui;
    }

    if (br.overrun()) {
        error = "truncated H.264 SPS";
        return false;
    }
    if (!s.flags.frame_mbs_only_flag && !s.flags.mb_adaptive_frame_field_flag) {
        // Pure field coding: every picture is a field, which we do not handle.
        error = "H.264 field-coded (interlaced) stream; only progressive and "
                "MBAFF frame pictures are supported";
        return false;
    }

    sps->max_frame_num = 1u << (s.log2_max_frame_num_minus4 + 4);
    sps->max_poc_lsb = 1u << (s.log2_max_pic_order_cnt_lsb_minus4 + 4);

    const uint32_t mbs_w = s.pic_width_in_mbs_minus1 + 1;
    const uint32_t mbs_h =
        (s.pic_height_in_map_units_minus1 + 1) * (s.flags.frame_mbs_only_flag ? 1 : 2);
    const int sub_w =
        (s.chroma_format_idc == STD_VIDEO_H264_CHROMA_FORMAT_IDC_420 ||
         s.chroma_format_idc == STD_VIDEO_H264_CHROMA_FORMAT_IDC_422) ? 2 : 1;
    const int sub_h = (s.chroma_format_idc == STD_VIDEO_H264_CHROMA_FORMAT_IDC_420) ? 2 : 1;
    const int crop_unit_x =
        (s.chroma_format_idc == STD_VIDEO_H264_CHROMA_FORMAT_IDC_MONOCHROME) ? 1 : sub_w;
    const int crop_unit_y =
        ((s.chroma_format_idc == STD_VIDEO_H264_CHROMA_FORMAT_IDC_MONOCHROME) ? 1 : sub_h) *
        (s.flags.frame_mbs_only_flag ? 1 : 2);
    sps->width = (int)(mbs_w * 16) -
                 crop_unit_x * (int)(s.frame_crop_left_offset + s.frame_crop_right_offset);
    sps->height = (int)(mbs_h * 16) -
                  crop_unit_y * (int)(s.frame_crop_top_offset + s.frame_crop_bottom_offset);

    const uint32_t frame_mbs = mbs_w * mbs_h;
    sps->max_dpb_frames =
        std::min<uint32_t>(std::max<uint32_t>(max_dpb_mbs(s.level_idc) / std::max(frame_mbs, 1u), 1), 16);
    sps->max_dpb_frames = std::max(sps->max_dpb_frames, (uint32_t)s.max_num_ref_frames);
    if (sps->max_num_reorder > sps->max_dpb_frames) sps->max_num_reorder = sps->max_dpb_frames;
    sps->valid = true;

    const bool first = !sps_[sps_id];
    sps_[sps_id] = std::move(sps);
    if (first || sps_id == 0) updateFormat(*sps_[sps_id]);
    have_sps_ = true;
    params_dirty_ = true;
    return true;
}

bool H264Decoder::parsePps(const uint8_t* nal, size_t size, std::string& error) {
    RbspReader br(nal, size, 1);
    auto pps = std::make_unique<Pps>();
    StdVideoH264PictureParameterSet& p = pps->std;

    const uint32_t pps_id = br.ue();
    const uint32_t sps_id = br.ue();
    if (pps_id >= 256 || sps_id >= 32) {
        error = "H.264 PPS/SPS id out of range";
        return false;
    }
    p.pic_parameter_set_id = (uint8_t)pps_id;
    p.seq_parameter_set_id = (uint8_t)sps_id;
    p.flags.entropy_coding_mode_flag = br.flag();
    p.flags.bottom_field_pic_order_in_frame_present_flag = br.flag();
    const uint32_t num_slice_groups = br.ue() + 1;
    if (num_slice_groups > 1) {
        error = "H.264 streams with slice groups (FMO) are not supported";
        return false;
    }
    p.num_ref_idx_l0_default_active_minus1 = (uint8_t)br.ue();
    p.num_ref_idx_l1_default_active_minus1 = (uint8_t)br.ue();
    p.flags.weighted_pred_flag = br.flag();
    p.weighted_bipred_idc = (StdVideoH264WeightedBipredIdc)br.u(2);
    p.pic_init_qp_minus26 = (int8_t)br.se();
    p.pic_init_qs_minus26 = (int8_t)br.se();
    p.chroma_qp_index_offset = (int8_t)br.se();
    p.flags.deblocking_filter_control_present_flag = br.flag();
    p.flags.constrained_intra_pred_flag = br.flag();
    p.flags.redundant_pic_cnt_present_flag = br.flag();

    if (br.moreRbspData()) {
        p.flags.transform_8x8_mode_flag = br.flag();
        p.flags.pic_scaling_matrix_present_flag = br.flag();
        if (p.flags.pic_scaling_matrix_present_flag) {
            const Sps* s = sps_[sps_id].get();
            const int chroma444 =
                (s && s->std.chroma_format_idc == STD_VIDEO_H264_CHROMA_FORMAT_IDC_444) ? 1 : 0;
            const int n_lists = 6 + (p.flags.transform_8x8_mode_flag ? (chroma444 ? 6 : 2) : 0);
            for (int i = 0; i < n_lists; ++i) {
                if (!br.flag()) continue;
                bool use_default = false;
                if (i < 6) {
                    parse_scaling_list(br, pps->scaling.ScalingList4x4[i], 16, use_default);
                } else {
                    parse_scaling_list(br, pps->scaling.ScalingList8x8[i - 6], 64, use_default);
                }
                pps->scaling.scaling_list_present_mask |= (uint16_t)(1u << i);
                if (use_default)
                    pps->scaling.use_default_scaling_matrix_mask |= (uint16_t)(1u << i);
            }
            p.pScalingLists = &pps->scaling;
        }
        p.second_chroma_qp_index_offset = (int8_t)br.se();
    } else {
        p.second_chroma_qp_index_offset = p.chroma_qp_index_offset;
    }
    if (br.overrun()) {
        error = "truncated H.264 PPS";
        return false;
    }
    pps->valid = true;
    pps_[pps_id] = std::move(pps);
    params_dirty_ = true;
    return true;
}

void H264Decoder::updateFormat(const Sps& sps) {
    format_.width = sps.width;
    format_.height = sps.height;
    format_.coded_width = (int)(sps.std.pic_width_in_mbs_minus1 + 1) * 16;
    format_.coded_height =
        (int)((sps.std.pic_height_in_map_units_minus1 + 1) *
              (sps.std.flags.frame_mbs_only_flag ? 1u : 2u)) * 16;
    format_.bit_depth = 8 + sps.std.bit_depth_luma_minus8;
    format_.chroma_format = (int)sps.std.chroma_format_idc;
    format_.max_dpb_slots = std::min<uint32_t>(sps.max_dpb_frames + 1, kMaxDpb);
    format_.max_active_references = std::min<uint32_t>(sps.max_dpb_frames, 16);
    format_.max_reorder = sps.max_num_reorder;
    format_.full_range = sps.vui.flags.video_full_range_flag != 0;
    format_.matrix_coefficients =
        sps.vui.flags.color_description_present_flag ? sps.vui.matrix_coefficients : 2;
    profile_.stdProfileIdc = sps.std.profile_idc;
    switch (sps.std.profile_idc) {
        case STD_VIDEO_H264_PROFILE_IDC_BASELINE: format_.profile_name = "Baseline"; break;
        case STD_VIDEO_H264_PROFILE_IDC_MAIN:     format_.profile_name = "Main"; break;
        case STD_VIDEO_H264_PROFILE_IDC_HIGH:     format_.profile_name = "High"; break;
        default:                                  format_.profile_name = "High 4:4:4"; break;
    }
}

// ---------------------------------------------------------------------------
// Slice header
// ---------------------------------------------------------------------------

bool H264Decoder::parseSliceHeader(const uint8_t* nal, size_t size, int nal_ref_idc, bool idr,
                                   SliceHeader& sh, std::string& error) {
    RbspReader br(nal, size, 1);
    sh = SliceHeader{};
    sh.first_mb = br.ue();
    sh.slice_type = br.ue() % 5;
    sh.pps_id = br.ue();
    if (sh.pps_id >= 256 || !pps_[sh.pps_id]) {
        error = "H.264 slice references unknown PPS " + std::to_string(sh.pps_id);
        return false;
    }
    const Pps& pps = *pps_[sh.pps_id];
    const uint32_t sps_id = pps.std.seq_parameter_set_id;
    if (!sps_[sps_id]) {
        error = "H.264 slice references unknown SPS " + std::to_string(sps_id);
        return false;
    }
    const Sps& sps = *sps_[sps_id];
    pending_sps_ = &sps;

    if (sps.std.flags.separate_colour_plane_flag) br.u(2);
    sh.frame_num = br.u(sps.std.log2_max_frame_num_minus4 + 4);
    if (!sps.std.flags.frame_mbs_only_flag) {
        sh.field_pic = br.flag();
        if (sh.field_pic) sh.bottom_field = br.flag();
    }
    if (sh.field_pic) {
        error = "H.264 field pictures are not supported (progressive/MBAFF only)";
        return false;
    }
    if (idr) sh.idr_pic_id = br.ue();
    if (sps.std.pic_order_cnt_type == STD_VIDEO_H264_POC_TYPE_0) {
        sh.poc_lsb = br.u(sps.std.log2_max_pic_order_cnt_lsb_minus4 + 4);
        if (pps.std.flags.bottom_field_pic_order_in_frame_present_flag && !sh.field_pic)
            sh.delta_poc_bottom = br.se();
    } else if (sps.std.pic_order_cnt_type == STD_VIDEO_H264_POC_TYPE_1 &&
               !sps.std.flags.delta_pic_order_always_zero_flag) {
        sh.delta_poc[0] = br.se();
        if (pps.std.flags.bottom_field_pic_order_in_frame_present_flag && !sh.field_pic)
            sh.delta_poc[1] = br.se();
    }
    if (pps.std.flags.redundant_pic_cnt_present_flag) br.ue();

    const bool is_b = sh.slice_type == 1;
    const bool is_p = sh.slice_type == 0 || sh.slice_type == 3;
    if (is_b) br.bit();  // direct_spatial_mv_pred_flag

    uint32_t num_l0 = pps.std.num_ref_idx_l0_default_active_minus1;
    uint32_t num_l1 = pps.std.num_ref_idx_l1_default_active_minus1;
    if (is_p || is_b) {
        if (br.flag()) {
            num_l0 = br.ue();
            if (is_b) num_l1 = br.ue();
        }
    }

    // ref_pic_list_modification(): skipped, but must be consumed to reach
    // dec_ref_pic_marking(). The driver re-parses it from the slice data.
    auto skip_modification = [&br]() {
        if (br.flag()) {
            while (true) {
                const uint32_t idc = br.ue();
                if (idc == 3 || br.overrun()) break;
                br.ue();
            }
        }
    };
    if (is_p || is_b) skip_modification();
    if (is_b) skip_modification();

    const int chroma_array =
        sps.std.flags.separate_colour_plane_flag ? 0 : (int)sps.std.chroma_format_idc;
    if ((pps.std.flags.weighted_pred_flag && is_p) ||
        (pps.std.weighted_bipred_idc == STD_VIDEO_H264_WEIGHTED_BIPRED_IDC_EXPLICIT && is_b)) {
        br.ue();  // luma_log2_weight_denom
        if (chroma_array != 0) br.ue();  // chroma_log2_weight_denom
        for (int list = 0; list < (is_b ? 2 : 1); ++list) {
            const uint32_t n = (list == 0 ? num_l0 : num_l1) + 1;
            for (uint32_t i = 0; i < n && !br.overrun(); ++i) {
                if (br.flag()) { br.se(); br.se(); }
                if (chroma_array != 0 && br.flag())
                    for (int j = 0; j < 2; ++j) { br.se(); br.se(); }
            }
        }
    }

    if (nal_ref_idc != 0) {
        if (idr) {
            sh.no_output_of_prior_pics = br.flag();
            sh.long_term_reference_flag = br.flag();
        } else {
            sh.adaptive_marking = br.flag();
            if (sh.adaptive_marking) {
                while (!br.overrun()) {
                    SliceHeader::Mmco m{};
                    m.op = br.ue();
                    if (m.op == 0) break;
                    if (m.op == 1 || m.op == 3) m.arg0 = br.ue();
                    if (m.op == 2) m.arg0 = br.ue();
                    if (m.op == 3 || m.op == 6) m.arg1 = br.ue();
                    if (m.op == 4) m.arg0 = br.ue();
                    sh.mmco.push_back(m);
                    if (sh.mmco.size() > 32) break;
                }
            }
        }
    }
    if (br.overrun()) {
        error = "truncated H.264 slice header";
        return false;
    }
    return true;
}

// ---------------------------------------------------------------------------
// Picture order count (8.2.1)
// ---------------------------------------------------------------------------

void H264Decoder::computePoc(const SliceHeader& sh, const Sps& sps, int nal_ref_idc, bool idr) {
    const int max_frame_num = (int)sps.max_frame_num;

    if (sps.std.pic_order_cnt_type == STD_VIDEO_H264_POC_TYPE_0) {
        int prev_msb = 0, prev_lsb = 0;
        if (!idr) {
            if (prev_mmco5_) {
                prev_msb = 0;
                prev_lsb = top_poc_;
            } else {
                prev_msb = prev_poc_msb_;
                prev_lsb = prev_poc_lsb_;
            }
        }
        const int max_lsb = (int)sps.max_poc_lsb;
        const int lsb = (int)sh.poc_lsb;
        int msb;
        if (lsb < prev_lsb && (prev_lsb - lsb) >= max_lsb / 2)
            msb = prev_msb + max_lsb;
        else if (lsb > prev_lsb && (lsb - prev_lsb) > max_lsb / 2)
            msb = prev_msb - max_lsb;
        else
            msb = prev_msb;
        top_poc_ = msb + lsb;
        bottom_poc_ = top_poc_ + sh.delta_poc_bottom;
        if (nal_ref_idc != 0) {
            prev_poc_msb_ = msb;
            prev_poc_lsb_ = lsb;
        }
    } else if (sps.std.pic_order_cnt_type == STD_VIDEO_H264_POC_TYPE_1) {
        int frame_num_offset;
        if (idr)
            frame_num_offset = 0;
        else if (prev_frame_num_ > (int)sh.frame_num)
            frame_num_offset = prev_frame_num_offset_ + max_frame_num;
        else
            frame_num_offset = prev_frame_num_offset_;

        const int cycle_len = sps.std.num_ref_frames_in_pic_order_cnt_cycle;
        int abs_frame_num = cycle_len != 0 ? frame_num_offset + (int)sh.frame_num : 0;
        if (nal_ref_idc == 0 && abs_frame_num > 0) --abs_frame_num;

        int expected = 0;
        if (abs_frame_num > 0 && cycle_len > 0) {
            int cycle_sum = 0;
            for (int i = 0; i < cycle_len; ++i) cycle_sum += sps.offset_for_ref_frame[i];
            const int cycle_cnt = (abs_frame_num - 1) / cycle_len;
            const int in_cycle = (abs_frame_num - 1) % cycle_len;
            expected = cycle_cnt * cycle_sum;
            for (int i = 0; i <= in_cycle; ++i) expected += sps.offset_for_ref_frame[i];
        }
        if (nal_ref_idc == 0) expected += sps.std.offset_for_non_ref_pic;
        top_poc_ = expected + sh.delta_poc[0];
        bottom_poc_ = top_poc_ + sps.std.offset_for_top_to_bottom_field + sh.delta_poc[1];
        prev_frame_num_offset_ = frame_num_offset;
    } else {
        int frame_num_offset;
        if (idr)
            frame_num_offset = 0;
        else if (prev_frame_num_ > (int)sh.frame_num)
            frame_num_offset = prev_frame_num_offset_ + max_frame_num;
        else
            frame_num_offset = prev_frame_num_offset_;
        int tmp;
        if (idr)
            tmp = 0;
        else if (nal_ref_idc == 0)
            tmp = 2 * (frame_num_offset + (int)sh.frame_num) - 1;
        else
            tmp = 2 * (frame_num_offset + (int)sh.frame_num);
        top_poc_ = bottom_poc_ = tmp;
        prev_frame_num_offset_ = frame_num_offset;
    }
    prev_frame_num_ = (int)sh.frame_num;
}

// ---------------------------------------------------------------------------
// Reference marking (8.2.5)
// ---------------------------------------------------------------------------

int H264Decoder::allocSlot() {
    const uint32_t n = std::max<uint32_t>(format_.max_dpb_slots, 1);
    for (uint32_t i = 0; i < n; ++i)
        if (!dpb_[i].used) return (int)i;
    // Every slot is a live reference: evict the oldest short-term one. A
    // conformant stream cannot get here, but a truncated file can.
    int worst = -1, worst_wrap = INT32_MAX;
    for (uint32_t i = 0; i < n; ++i)
        if (!dpb_[i].long_term && dpb_[i].frame_num_wrap < worst_wrap) {
            worst_wrap = dpb_[i].frame_num_wrap;
            worst = (int)i;
        }
    if (worst < 0) worst = 0;
    dpb_[worst].used = false;
    return worst;
}

void H264Decoder::applyMarking(const SliceHeader& sh, const Sps& sps) {
    const int max_frame_num = (int)sps.max_frame_num;
    const int curr_pic_num = (int)sh.frame_num;

    // Refresh FrameNumWrap / PicNum for the existing short-term references.
    for (auto& f : dpb_) {
        if (!f.used || f.long_term) continue;
        f.frame_num_wrap = (f.frame_num > curr_pic_num) ? f.frame_num - max_frame_num
                                                         : f.frame_num;
        f.pic_num = f.frame_num_wrap;
    }

    if (pending_idr_) {
        for (auto& f : dpb_) f.used = false;
        return;
    }
    if (!pending_is_ref_) return;

    if (sh.adaptive_marking) {
        for (const auto& m : sh.mmco) {
            switch (m.op) {
                case 1: {
                    const int pic_num_x = curr_pic_num - (int)(m.arg0 + 1);
                    for (auto& f : dpb_)
                        if (f.used && !f.long_term && f.pic_num == pic_num_x) f.used = false;
                    break;
                }
                case 2:
                    for (auto& f : dpb_)
                        if (f.used && f.long_term && f.long_term_frame_idx == (int)m.arg0)
                            f.used = false;
                    break;
                case 3: {
                    const int pic_num_x = curr_pic_num - (int)(m.arg0 + 1);
                    for (auto& f : dpb_)
                        if (f.used && f.long_term && f.long_term_frame_idx == (int)m.arg1)
                            f.used = false;
                    for (auto& f : dpb_)
                        if (f.used && !f.long_term && f.pic_num == pic_num_x) {
                            f.long_term = true;
                            f.long_term_frame_idx = (int)m.arg1;
                            f.std_ref.flags.used_for_long_term_reference = 1;
                        }
                    break;
                }
                case 4: {
                    const int max_idx = (int)m.arg0 - 1;
                    for (auto& f : dpb_)
                        if (f.used && f.long_term && f.long_term_frame_idx > max_idx)
                            f.used = false;
                    break;
                }
                case 5:
                    for (auto& f : dpb_) f.used = false;
                    prev_mmco5_ = true;
                    break;
                case 6:
                    for (auto& f : dpb_)
                        if (f.used && f.long_term && f.long_term_frame_idx == (int)m.arg1)
                            f.used = false;
                    break;
                default: break;
            }
        }
        return;
    }

    // Sliding window: drop the oldest short-term reference once the buffer is
    // at max_num_ref_frames.
    int n_short = 0, n_long = 0;
    for (const auto& f : dpb_) {
        if (!f.used) continue;
        (f.long_term ? n_long : n_short)++;
    }
    const int max_refs = std::max<int>(sps.std.max_num_ref_frames, 1);
    if (n_short + n_long >= max_refs) {
        int oldest = -1, oldest_wrap = INT32_MAX;
        for (uint32_t i = 0; i < kMaxDpb; ++i)
            if (dpb_[i].used && !dpb_[i].long_term && dpb_[i].frame_num_wrap < oldest_wrap) {
                oldest_wrap = dpb_[i].frame_num_wrap;
                oldest = (int)i;
            }
        if (oldest >= 0) dpb_[oldest].used = false;
    }
}

// ---------------------------------------------------------------------------
// CodecDecoder
// ---------------------------------------------------------------------------

bool H264Decoder::init(const TrackInfo& track, std::string& error) {
    format_.width = track.width;
    format_.height = track.height;
    const std::vector<uint8_t>& cfg = track.codec_config;
    if (cfg.size() >= 7) {
        // AVCDecoderConfigurationRecord.
        size_t p = 5;
        const int n_sps = cfg[p++] & 0x1f;
        for (int i = 0; i < n_sps && p + 2 <= cfg.size(); ++i) {
            const size_t len = ((size_t)cfg[p] << 8) | cfg[p + 1];
            p += 2;
            if (p + len > cfg.size()) break;
            if (!parseSps(cfg.data() + p, len, error)) return false;
            p += len;
        }
        if (p < cfg.size()) {
            const int n_pps = cfg[p++];
            for (int i = 0; i < n_pps && p + 2 <= cfg.size(); ++i) {
                const size_t len = ((size_t)cfg[p] << 8) | cfg[p + 1];
                p += 2;
                if (p + len > cfg.size()) break;
                if (!parsePps(cfg.data() + p, len, error)) return false;
                p += len;
            }
        }
    }
    if (!have_sps_) {
        error = "H.264 track carries no SPS in its configuration record; in-band "
                "parameter sets are supported but the session cannot be sized "
                "without one";
        return false;
    }
    return true;
}

VkVideoSessionParametersKHR H264Decoder::createParameters(VkDevice device,
                                                           VkVideoSessionKHR session,
                                                           std::string& error) {
    std::vector<StdVideoH264SequenceParameterSet> sps_list;
    std::vector<StdVideoH264PictureParameterSet>  pps_list;
    for (const auto& s : sps_)
        if (s && s->valid) sps_list.push_back(s->std);
    for (const auto& p : pps_)
        if (p && p->valid) pps_list.push_back(p->std);

    VkVideoDecodeH264SessionParametersAddInfoKHR add{
        VK_STRUCTURE_TYPE_VIDEO_DECODE_H264_SESSION_PARAMETERS_ADD_INFO_KHR};
    add.stdSPSCount = (uint32_t)sps_list.size();
    add.pStdSPSs = sps_list.data();
    add.stdPPSCount = (uint32_t)pps_list.size();
    add.pStdPPSs = pps_list.data();

    VkVideoDecodeH264SessionParametersCreateInfoKHR codec{
        VK_STRUCTURE_TYPE_VIDEO_DECODE_H264_SESSION_PARAMETERS_CREATE_INFO_KHR};
    codec.maxStdSPSCount = 32;
    codec.maxStdPPSCount = 256;
    codec.pParametersAddInfo = &add;

    VkVideoSessionParametersCreateInfoKHR ci{
        VK_STRUCTURE_TYPE_VIDEO_SESSION_PARAMETERS_CREATE_INFO_KHR};
    ci.pNext = &codec;
    ci.videoSession = session;

    VkVideoSessionParametersKHR out = VK_NULL_HANDLE;
    const VkResult r = video_api().createParameters(device, &ci, nullptr, &out);
    if (r != VK_SUCCESS) {
        error = "vkCreateVideoSessionParametersKHR (H.264) failed";
        return VK_NULL_HANDLE;
    }
    params_dirty_ = false;
    return out;
}

void H264Decoder::flush() {
    for (auto& f : dpb_) f.used = false;
    prev_poc_msb_ = prev_poc_lsb_ = 0;
    prev_frame_num_ = prev_frame_num_offset_ = 0;
    prev_mmco5_ = false;
    pending_slot_ = -1;
}

bool H264Decoder::decodeFrame(const uint8_t* data, size_t size, int nal_length_size,
                              std::vector<uint8_t>& bitstream,
                              std::vector<uint32_t>& slice_offsets, PictureInfo& out,
                              std::string& error) {
    bitstream.clear();
    slice_offsets.clear();
    out = PictureInfo{};

    bool have_slice = false;
    int  nal_ref_idc = 0;
    bool idr = false;
    SliceHeader sh{};

    // Walk the access unit's NAL units. `nal_length_size == 0` means the data
    // is already Annex-B.
    size_t p = 0;
    while (p < size) {
        const uint8_t* nal = nullptr;
        size_t nal_size = 0;
        if (nal_length_size > 0) {
            if (p + (size_t)nal_length_size > size) break;
            size_t len = 0;
            for (int i = 0; i < nal_length_size; ++i) len = (len << 8) | data[p + i];
            p += (size_t)nal_length_size;
            if (len == 0 || p + len > size) break;
            nal = data + p;
            nal_size = len;
            p += len;
        } else {
            // Annex-B: find this start code and the next.
            while (p + 3 <= size && !(data[p] == 0 && data[p + 1] == 0 &&
                                      (data[p + 2] == 1 ||
                                       (p + 4 <= size && data[p + 2] == 0 && data[p + 3] == 1))))
                ++p;
            if (p + 3 > size) break;
            size_t start = (data[p + 2] == 1) ? p + 3 : p + 4;
            size_t q = start;
            while (q + 3 <= size &&
                   !(data[q] == 0 && data[q + 1] == 0 && (data[q + 2] == 1 || data[q + 2] == 0)))
                ++q;
            if (q + 3 > size) q = size;
            nal = data + start;
            nal_size = q - start;
            p = q;
        }
        if (nal_size == 0) continue;

        const int type = nal[0] & 0x1f;
        const int ref_idc = (nal[0] >> 5) & 3;
        switch (type) {
            case 7:
                if (!parseSps(nal, nal_size, error)) return false;
                break;
            case 8:
                if (!parsePps(nal, nal_size, error)) return false;
                break;
            case 1:
            case 5: {
                if (!have_slice) {
                    nal_ref_idc = ref_idc;
                    idr = (type == 5);
                    if (!parseSliceHeader(nal, nal_size, ref_idc, idr, sh, error)) return false;
                    have_slice = true;
                }
                slice_offsets.push_back((uint32_t)bitstream.size());
                annexb_append(bitstream, nal, nal_size);
                break;
            }
            default:
                break;  // SEI, AUD, filler: the driver does not need them
        }
    }

    if (!have_slice) {
        error = "H.264 access unit contains no coded slice";
        return false;
    }
    const Sps& sps = *pending_sps_;

    prev_mmco5_ = false;
    computePoc(sh, sps, nal_ref_idc, idr);
    if (idr) {
        for (auto& f : dpb_) f.used = false;
    }

    pending_sh_ = sh;
    pending_is_ref_ = nal_ref_idc != 0;
    pending_idr_ = idr;

    // Reference slots active for this picture, and the slot it will occupy.
    // The marking process runs in commitFrame(), so the slot list here is the
    // state as of the previous picture, which is exactly what the driver needs.
    out.refs.clear();
    for (uint32_t i = 0; i < format_.max_dpb_slots; ++i) {
        if (!dpb_[i].used) continue;
        dpb_[i].std_ref.FrameNum = (uint16_t)dpb_[i].frame_num;
        dpb_[i].std_ref.PicOrderCnt[0] = dpb_[i].poc;
        dpb_[i].std_ref.PicOrderCnt[1] = dpb_[i].poc;
        out.refs.push_back({(int32_t)i, &dpb_[i].std_ref});
        if (out.refs.size() >= format_.max_active_references) break;
    }

    const int slot = allocSlot();
    pending_slot_ = slot;

    std_pic_ = StdVideoDecodeH264PictureInfo{};
    std_pic_.seq_parameter_set_id = sps.std.seq_parameter_set_id;
    std_pic_.pic_parameter_set_id = (uint8_t)sh.pps_id;
    std_pic_.frame_num = (uint16_t)sh.frame_num;
    std_pic_.idr_pic_id = (uint16_t)sh.idr_pic_id;
    std_pic_.PicOrderCnt[0] = top_poc_;
    std_pic_.PicOrderCnt[1] = bottom_poc_;
    std_pic_.flags.field_pic_flag = 0;
    std_pic_.flags.is_intra = (sh.slice_type == 2 || sh.slice_type == 4) ? 1 : 0;
    std_pic_.flags.IdrPicFlag = idr ? 1 : 0;
    std_pic_.flags.bottom_field_flag = 0;
    std_pic_.flags.is_reference = pending_is_ref_ ? 1 : 0;
    std_pic_.flags.complementary_field_pair = 0;

    vk_pic_ = VkVideoDecodeH264PictureInfoKHR{
        VK_STRUCTURE_TYPE_VIDEO_DECODE_H264_PICTURE_INFO_KHR};
    vk_pic_.pStdPictureInfo = &std_pic_;
    vk_pic_.sliceCount = (uint32_t)slice_offsets.size();
    vk_pic_.pSliceOffsets = slice_offsets.data();

    setup_ref_ = StdVideoDecodeH264ReferenceInfo{};
    setup_ref_.FrameNum = (uint16_t)sh.frame_num;
    setup_ref_.PicOrderCnt[0] = top_poc_;
    setup_ref_.PicOrderCnt[1] = bottom_poc_;
    setup_ref_.flags.used_for_long_term_reference =
        (idr && sh.long_term_reference_flag) ? 1 : 0;

    out.decode_pnext = &vk_pic_;
    // Every picture activates a slot: the driver writes the reconstructed frame
    // there, and a non-reference picture simply has its slot freed again in
    // commitFrame().
    out.setup_slot = slot;
    out.setup_std_ref = &setup_ref_;
    out.poc = std::min(top_poc_, bottom_poc_);
    out.params_changed = params_dirty_;
    out.sequence_restart = idr;
    return true;
}

void H264Decoder::commitFrame() {
    if (pending_slot_ < 0 || !pending_sps_) return;
    applyMarking(pending_sh_, *pending_sps_);

    DpbFrame& f = dpb_[(size_t)pending_slot_];
    f = DpbFrame{};
    if (pending_is_ref_) {
        f.used = true;
        f.frame_num = (int)pending_sh_.frame_num;
        f.frame_num_wrap = f.frame_num;
        f.pic_num = f.frame_num;
        f.poc = std::min(top_poc_, bottom_poc_);
        f.long_term = pending_idr_ && pending_sh_.long_term_reference_flag;
        if (!pending_idr_) {
            for (const auto& m : pending_sh_.mmco)
                if (m.op == 6) {
                    f.long_term = true;
                    f.long_term_frame_idx = (int)m.arg1;
                }
        }
        f.std_ref.flags.used_for_long_term_reference = f.long_term ? 1 : 0;
        f.std_ref.FrameNum = (uint16_t)f.frame_num;
        f.std_ref.PicOrderCnt[0] = f.poc;
        f.std_ref.PicOrderCnt[1] = f.poc;
    }
    if (prev_mmco5_) {
        // 8.2.1: after MMCO 5 the picture's POC is rebased to zero.
        bottom_poc_ -= top_poc_;
        top_poc_ = 0;
        prev_frame_num_ = 0;
        prev_frame_num_offset_ = 0;
    }
    pending_slot_ = -1;
}

}  // namespace

// Three-byte start codes, not four: the offsets handed to the driver in
// pSliceOffsets / pSliceSegmentOffsets point at the start code, and NVIDIA's
// H.265 slice-header parser assumes a 3-byte prefix there -- a 4-byte one
// leaves it reading the slice header one byte early.
void annexb_append(std::vector<uint8_t>& out, const uint8_t* nal, size_t size) {
    out.push_back(0);
    out.push_back(0);
    out.push_back(1);
    out.insert(out.end(), nal, nal + size);
}

std::unique_ptr<CodecDecoder> make_h264_decoder() { return std::make_unique<H264Decoder>(); }

}  // namespace video
