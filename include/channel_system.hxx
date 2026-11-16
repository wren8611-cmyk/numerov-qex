// channel_system.hxx

#pragma once

#include "channel_info.hxx"

#include <iostream>
#include <vector>
#include <unordered_map>
#include <stdexcept>
#include <cstddef>
#include <utility>

class ChannelSystem {
public:
    // Constructor
    explicit ChannelSystem(std::vector<ChannelInfo> channels) : channels_(std::move(channels))
    {
        if (channels_.empty()) throw std::invalid_argument("Channel system cannot be empty");

        for (std::size_t i = 0; i < channels_.size(); ++i) {
            if (!channel_map_.emplace(channels_[i].id, i).second) throw std::invalid_argument("Duplicate channel ID");
        }
    }

    // Accessors
    std::size_t nchannels() const { return channels_.size(); }

    const std::vector<ChannelInfo>& channels()   const { return channels_; }
    const ChannelInfo& operator[](std::size_t i) const { return channels_.at(i); }

    std::size_t channel_index(int id) const { return channel_map_.at(id); }
    const ChannelInfo& channel_by_id(int id) const { return channels_.at(channel_index(id)); }

    // Display channel system information
    void info(std::ostream& os = std::cout) const {
        os << "Channel System Information\n";
        os << "==========================\n";
        os << "Number of channels : " << nchannels() << "\n\n";

        for (std::size_t i = 0; i < nchannels(); ++i) {
            os << "Matrix index : " << i << '\n';
            channels_[i].summary(os);
            if (i + 1 < nchannels()) os << '\n';
        }
    }

private:
    std::vector<ChannelInfo> channels_;
    std::unordered_map<int, std::size_t> channel_map_;
};