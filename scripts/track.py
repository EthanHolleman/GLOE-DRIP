import pandas as pd
from pathlib import Path


class FilepathToTrackConverter():
    
    def __init__(self, filepath, url_prefix=''):
        self.filepath = Path(filepath)
        self.bigDataUrl = url_prefix + self.filepath.name
  

class gloeRepConverter(FilepathToTrackConverter):
    
    def __init__(self, filepath, url_prefix):
        super().__init__(filepath, url_prefix)
    
    @property
    def color(self):
        if 'fwd' in self.filepath.name:
            return '66, 135, 245'
        else:
            return '245, 66, 66'
    
    @property
    def shortLabel(self):
        return self.filepath.name.replace('_', ' ')
    
    @property
    def sample_name(self):
        return self.filepath.name.split('.')[0]
    
    @property
    def modification(self):
        return str(
            snakemake.params.samples.loc[self.sample_name]['modification']
            ).replace(' ', '_')
    
    
    def make_track_string(self):
        ts = f'''
    track {self.filepath.stem}
    bigDataUrl {self.bigDataUrl}
    shortLabel {self.shortLabel}
    longLabel {self.shortLabel}
    type bigWig
    parent {self.modification}
    color {self.color}
    autoScale off
    maxHeightPixels 40:40:40
    viewLimits 0.25:7
    windowingFunction mean
    compositeTrack on
    subGroup1 view View Signal=Signal Peaks=Peaks
    visibility full
    trackDescription {self.filepath.stem}.html
    priority 
    subGroups view=Signal\n
        '''
        return ts


class gloePeakConverter(gloeRepConverter):

    def __init__(self, filepath, url_prefix):
        super().__init__(filepath, url_prefix)
    
    @property
    def parent(self):
        return self.filepath.stem.split('.')[0]

    
    def make_track_string(self):
        ts = f'''
    track {self.filepath.stem}
    bigDataUrl {self.bigDataUrl}
    shortLabel {self.shortLabel}
    longLabel {self.shortLabel}
    type bigWig
    color {self.color}
    autoScale off
    maxHeightPixels 40:40:40
    viewLimits 0.25:7
    windowingFunction mean
    compositeTrack on
    subGroup1 view View Signal=Signal Peaks=Peaks
    visibility full
    trackDescription {self.filepath.stem}.html
    priority 1
    subGroups view=Signal\n
        '''
        return ts
    


def add_parent_tracks(samples):
    track_strings = []
    for row_index, row in samples.iterrows():
        ts = f'''
    track {str(row['modification']).replace(' ', '_')}
    type bigWig
    view Signal
    visibility full
    shortLabel {str(row['modification']).replace(' ', '_')}
    longLabel {str(row['modification']).replace(' ', '_')}
    autoScale off
    maxHeightPixels 40:40:40
    viewLimits 0.25:7
    windowingFunction mean
    compositeTrack on
    subGroup1 view View Signal=Signal Peaks=Peaks
    visibility full
    trackDescription {str(row['modification']).replace(' ', '_')}.html
    priority 1
    '''
    #     track_strings.append(ts)
    # for call_type in ['CAT', 'TAT']:
    #     ts = f'''
    # track {call_type}
    # type bigWig
    # view Signal
    # visibility full
    # shortLabel {call_type}
    # longLabel {call_type}
    # autoScale off
    # maxHeightPixels 40:40:40
    # viewLimits 0.25:7
    # windowingFunction mean
    # compositeTrack on
    # subGroup1 view View Signal=Signal Peaks=Peaks
    # visibility full
    # trackDescription {call_type}.html
    # priority 1
    # '''
        track_strings.append(ts)
    #parent = ['track GLOE-seq_replicates\ncompositeTrack on\nsubGroup1 view View Signal=Signal Peaks=Peaks\nvisibility full\npriority 1\nshortLabel GLOE-seq Replicates\nlongLabel GLOE-seq Replicates\ntype bigWig\ndescriptionUrl GLOE-seq_replicates.html\n']
    return list(set(track_strings))


def write_tracks_to_output(track_objects, output_file):
    with open(str(output_file), 'a') as handle:
        for track in track_objects:
            handle.write(track.make_track_string())
    return output_file


def write_parent_tracks(parent_track_strings, output_file):
    with open(str(output_file), 'w') as handle:
        for track in parent_track_strings:
            handle.write(track)
    return output_file


def input_files_to_tracks():
    # assume snakemake context
    tracks = []
    for input_file in snakemake.input.reps:
        tracks.append(gloeRepConverter(input_file, snakemake.params.url))
    for input_file in snakemake.input.peaks:
        print(input_file)
        tracks.append(gloePeakConverter(input_file, snakemake.params.url))
    return tracks


def main():
    print('starting')
    print(snakemake.input.peaks)
    print('peaks')
    tracks = input_files_to_tracks()
    parent_tracks = add_parent_tracks(snakemake.params.samples)
    write_parent_tracks(parent_tracks, snakemake.output)
    write_tracks_to_output(tracks, snakemake.output)

if __name__ == '__main__':
    main()

    




    
    
    
    
    

    



# class track():


#     def __init__(self, track, track_data=None, child_tracks=[],
#                 priority=None, shortLabel=None, longLabel=None,
#                 view='Signal', visibility='full', autoscale='off')
#                 **kwargs):
#         self.track = track  # name of track
#         self.child_tracks = child_tracks
#         self.parent = parent
#         # update with kwargs
    
        
#     def make_track_string(self, indent=0):
#         attributes_sorted = sorted(list(self.__dict__.keys()))
#         indent_str = '\t' * indent
#         attribute_string = '\n'.join(
#             [f'{indent_str}{key} {self.__dict__[key]}' for key in attributes_sorted])
        

    
#     def complete_track_file(self, output_path):  # assume have parent track
#         # get the head track
#         with open(output_path, 'w') as handle:
#             def add_track_to_trackfile(tracks):
#                 for track in tracks:
#                     handle.write(track.make_track_string())
#                     if track.child_tracks:
#                         add_track_to_trackfile(track.child_tracks)
#         return output_path
        
        


            
        


        
        



