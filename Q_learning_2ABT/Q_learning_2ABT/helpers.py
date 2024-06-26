import numpy as np
import pickle
from pathlib import Path
import pandas as pd
import torch
import matplotlib.pyplot as plt


def event_triggered_traces(arr, trigger_signal, win_bounds, trigger_signal_is_idx=False):
    '''
    Makes event triggered traces along last dimension
    RH 2021
    
    Args:
        arr (np.ndarray):
            Input array. Last dimension will be 
             aligned to boolean True values in 
             'trigger_signal'
        trigger_signal (boolean np.ndarray):
            1-D boolean array. True values are trigger
             events
        win_bounds (size 2 integer list or np.ndarray):
            2 value integer array. win_bounds[0] should
             be negative and is the number of samples prior
             to the event that the window starts. 
             win_bounds[1] is the number of samples 
             following the event.
            Events that would have a window extending
             before or after the bounds of the length
             of the trace are discarded.
        trigger_signal_is_idx (bool):
            If True then trigger_signal is an index array.
            If False then trigger_signal is a boolean array.
            Use an index array if there are multiple events
             with the same index, else they will be
             collapsed when this is 'True'.
     
     Returns:
        et_traces (np.ndarray):
             Event Triggered Traces. et_traces.ndim = 
              arr.ndim+1. Last dimension is new and is
              the event number axis. Note that events 
              near edges are discarded if window extends
              past edge bounds.
        xAxis (np.ndarray):
            x-axis of the traces. Aligns with dimension
             et_traces.shape[-2]
        windows (np.ndarray):
    '''
    def bounds_to_win(x_pos, win_bounds):
        return x_pos + np.arange(win_bounds[0], win_bounds[1])
    def make_windows(x_pos, win_bounds):
        return np.apply_along_axis(bounds_to_win, 0, tuple([x_pos]), win_bounds).T

    axis=arr.ndim-1
    len_axis = arr.shape[axis]

    if trigger_signal_is_idx:
        windows = make_windows(trigger_signal[np.isnan(trigger_signal)==0], win_bounds)
    else:
        windows = make_windows(np.nonzero(trigger_signal)[0], win_bounds)
    win_toInclude = (np.sum(windows<0, axis=1)==0) * (np.sum(windows>len_axis, axis=1)==0)
    win_toExclude = win_toInclude==False
    n_win_included = np.sum(win_toInclude)
    n_win_excluded = np.sum(win_toExclude)
    windows = windows[win_toInclude]


    windows_flat = np.reshape(windows, (windows.size))
    
    axes_all = np.arange(arr.ndim)
    axes_non = axes_all[axes_all != axis]
    et_traces_flat = np.take_along_axis(arr, np.expand_dims(np.array(windows_flat, dtype=np.int64), axis=tuple(axes_non)), axis=axis)

    new_shape = np.array(et_traces_flat.shape)
    new_shape[axis] = new_shape[axis] / windows.shape[1]
    new_shape = np.concatenate((new_shape, np.array([windows.shape[1]])))
    et_traces = np.reshape(et_traces_flat, new_shape)
    
    xAxis = np.arange(win_bounds[0], win_bounds[1])
    
    return et_traces, xAxis, windows

#####################################################
# saving and loading pickles
def pickle_save(obj, path_save, mode='wb', mkdir=False, allow_overwrite=True):
    Path(path_save).parent.mkdir(parents=True, exist_ok=True) if mkdir else None
    assert allow_overwrite or not Path(path_save).exists(), f'{path_save} already exists.'
    assert Path(path_save).parent.exists(), f'{Path(path_save).parent} does not exist.'
    assert Path(path_save).parent.is_dir(), f'{Path(path_save).parent} is not a directory.'
    with open(path_save, mode) as f:
        pickle.dump(obj, f,)

def pickle_load(filename, mode='rb'):
    with open(filename, mode) as f:
        return pickle.load(f)
#######################################################
    

# add a key value pair to a dictionary       
def append_dict(d, key, val):
    d2 = d.copy()
    d2[key] = val
    return d2

# turn our logger output into a dataframe
def logger_to_df(l):
    l2 = {key: val.detach().cpu().numpy() if isinstance(val, torch.Tensor) else val for key, val in l.items()}
    return pd.DataFrame(l2)

########################################################
def idx_to_oneHot(arr, n_classes=None, dtype=None):
    """
    Convert an array of class indices to matrix of
     one-hot vectors.
    RH 2021
    Args:
        arr (np.ndarray):
            1-D array of class indices.
            Values should be integers >= 0.
            Values will be used as indices in the
             output array.
        n_classes (int):
            Number of classes.
    
    Returns:
        oneHot (np.ndarray):
            2-D array of one-hot vectors.
    """
    if type(arr) is np.ndarray:
        max = np.max
        zeros = np.zeros
        arange = np.arange
        dtype = np.bool8 if dtype is None else dtype
    elif type(arr) is torch.Tensor:
        max = torch.max
        zeros = torch.zeros
        arange = torch.arange
        dtype = torch.bool if dtype is None else dtype
    assert arr.ndim == 1

    if n_classes is None:
        n_classes = max(arr)+1
    oneHot = zeros((len(arr), n_classes), dtype=dtype)
    oneHot[arange(len(arr)), arr] = True
    return oneHot

def confusion_matrix(y_hat, y_true):
    """
    Compute the confusion matrix from y_hat and y_true.
    y_hat should be either predictions or probabilities.
    RH 2021
    Args:
        y_hat (np.ndarray): 
            numpy array of predictions or probabilities. 
            Either PREDICTIONS: 2-D array of booleans
             ('one hots') or 1-D array of predicted 
             class indices.
            Or PROBABILITIES: 2-D array floats ('one hot
             like')
        y_true (np.ndarray):
            Either 1-D array of true class indices OR a
             precomputed onehot matrix.
    Returns:
        conf_mat (np.ndarray):
            2-D array of confusion matrix.
            Columns are true classes, rows are predicted.
            Note that this is the transpose of the
             sklearn convention.
    """
    n_classes = max(np.max(y_true)+1, np.max(y_hat)+1)
    if y_hat.ndim == 1:
        y_hat = idx_to_oneHot(y_hat, n_classes)
    if y_true.ndim == 1:
        y_true = idx_to_oneHot(y_true, n_classes)
    cmat = np.dot(y_hat.T, y_true)
    return cmat / np.sum(cmat, axis=0)[None,:]



def save_fig_types(fig, base_folder, fn, types=['png', 'svg', 'pdf', 'eps']):
    """
    Saves figures as multiple different types.
    fig: matplotlib.pyplot.figure
        matplotlib figure to save
    base_fn: str
        base_filename (including directory) to which to append extension
    types:
        list of file extensions to include as figure saves
    """
    base_path = Path(base_folder)
    base_path.mkdir(parents=True, exist_ok=True)
    for typ in types:
        base_file = str((base_path / f'{fn}.{typ}').resolve())
        fig.savefig(base_file)
    return