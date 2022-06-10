from sklearn import linear_model
# Importing OrdinalEncoder from Sklearn
# library from preprocessing Module
from sklearn.preprocessing import OrdinalEncoder


def make_design_matrix(sbj, drop_log, task_sep=True):
    '''

        This function generates a design matrix that corresponds to features of each trial in a groove task.
        The features included are intercepts, rhythm, harmony, and their interactions, rating, task, and intercept.

        Parameters:
            sbj (String) : Subject Identifier
            drop_log (MNE object): List of tuples that is the same length as the number of epochs. 
                                   Each tuple indicates whether the epoch should be dropped for analysis and why.

        Returns:
            design_matrix (ndarray) : Array of shape (n_trials, n_features). 
            features (list) : Names of features.
    '''

    # Creating a instance of label Encoder.
    # We only need one for harmony and rhythm
    # because harmony set = {'high', 'medium'}
    # and rhythm set = {'high', 'medium'} 
    harm_rhyth_le = OrdinalEncoder()
    
    # Structure of behavioral data
    # First two fields are cell arrays that have rhythm, and harmony
    # Third field is the choice for that trial
    # fourth field is the response time for that trial
    # Fifth field is the trl(triggers)

    for block in range(1,3):
        behavioral_data_path = f"/home/knight/groove/data/{sbj}/{sbj}_behavior_%03d_block.mat" % block
        behavioral_data = scipy.io.loadmat(behavioral_data_path)

        rhythms = np.unique(behavioral_data['behavioraldata'][0,0]['rhythm'])
        harmonies = np.unique(behavioral_data['behavioraldata'][0,0]['harmony'])
#         choices = np.unique(behavioral_data['behavioraldata'][0,0]['choices'])

        if block == 1:           
            offset = 0
            future_offset = len(behavioral_data['behavioraldata'][0,0]['choices'])
        else:
            offset = future_offset-1
            future_offset += len(behavioral_data['behavioraldata'][0,0]['choices'])

        # We drop all information about epochs that were dropped
        filtered_rhythms = [rhythm[0][0] for i, rhythm in 
                            enumerate(behavioral_data['behavioraldata'][0,0]['rhythm']) if not drop_log[i+offset]]
        filtered_harmonies = [harmony[0][0] for i, harmony in 
                            enumerate(behavioral_data['behavioraldata'][0,0]['harmony']) if not drop_log[i+offset]]
        filtered_choices = [choice[0] for i, choice in 
                            enumerate(behavioral_data['behavioraldata'][0,0]['choices']) if not drop_log[i+offset]]
        print(f"choice {len(behavioral_data['behavioraldata'][0,0]['choices'])} ")

        center_ratings = 3
        # For groove experiment, block corresponds to different task
        if block == 1:
            # possible_single_sorts = ['Choices', 'Rhythms', 'Harmonies']
            # sorts_ind = [choices_sorted_indices, rhythms_sorted_indices, harmonies_sorted_indices]
            # filtered_things = [filtered_choices, filtered_rhythms, filtered_harmonies]
            harmony_labels = harm_rhyth_le.fit_transform(np.array([filtered_harmonies]).T)
            if task_sep:
                pleasure_ratings = np.zeros((len(filtered_choices),1))
                groove_ratings = np.array(filtered_choices).reshape(len(filtered_choices),1).astype('int')-center_ratings
            else:
                task = np.ones((len(filtered_choices),1))-1
        else:
            harmony_labels = harm_rhyth_le.transform(np.array([filtered_harmonies]).T)
            if task_sep:
                pleasure_ratings = np.array(filtered_choices).reshape(len(filtered_choices),1).astype('int')-center_ratings
                groove_ratings = np.zeros((len(filtered_choices),1))
            else:
                task = np.ones((len(filtered_choices),1))

        # For groove experiment, harmony is medium and high, and rhythm is also medium and high so we can use the same encoder
        rhythm_labels = harm_rhyth_le.transform(np.array([filtered_rhythms]).T)

        # Encoder makes binary, but we want -1,1 binary not 0,1 binary
        harmony_labels[harmony_labels <= 0] = -1.
        rhythm_labels[rhythm_labels <= 0] = -1.

        # Interaction term is simpler because of -1,1 binary encoding
        interaction_term = harmony_labels*rhythm_labels


        ratings = np.array(filtered_choices).reshape(len(filtered_choices),1)
        if task_sep:
            features=['Harmony','Rhythm','Interaction Term', 'Groove', 'Pleasure']
        else:
            features = ['Harmony', 'Rhythm', 'Interaction Term', 'Ratings', 'Task']

        # Important! All ratings are centered at 3 so mean ~ 0 (not really because more patients don't use the full scale)
        design_matrix_block = np.hstack((harmony_labels, rhythm_labels, interaction_term, ratings.astype('int')-center_ratings, task))

        if block == 1:
            design_matrix = design_matrix_block
        else:
            design_matrix = np.vstack((design_matrix, design_matrix_block))
        print(design_matrix.shape)
    return design_matrix, features