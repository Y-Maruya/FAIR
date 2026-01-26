#pragma once


struct TrackFindAlgOutput
{
    bool found = false; // Whether a track was found
    int num_tracks = 0; // Number of tracks found
    // settings for track following
    int isolation_distance = 42;//mm per layer (20mm)
    int min_hits_in_track = 6;
    int max_missing_layers = 2;
    int max_distance_xy = 42;//mm
    int num_later_layers_to_fit = 10;
};

inline std::vector<FieldDesc> describe(const SimpleFittedTrack*) {
    return {
        field("found", &TrackFindAlgOutput::found),
        field("num_tracks", &TrackFindAlgOutput::num_tracks),
        field("isolation_distance", &TrackFindAlgOutput::isolation_distance),
        field("min_hits_in_track", &TrackFindAlgOutput::min_hits_in_track),
        field("max_missing_layers", &TrackFindAlgOutput::max_missing_layers),
        field("max_distance_xy", &TrackFindAlgOutput::max_distance_xy),
        field("num_later_layers_to_fit", &TrackFindAlgOutput::num_later_layers_to_fit),
    };
}

