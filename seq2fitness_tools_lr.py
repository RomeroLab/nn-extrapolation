# Import necessary libraries here
from os.path import abspath
import sys

module_path = abspath("nn4dms_nn-extrapolate/code")
if module_path not in sys.path:
    sys.path.append(module_path)

print(sys.path)

import numpy as np
import encoding
import inference_lr as inference  # use this inference for LR model only - removed training_ph
# import inference  # use this inference for any model that is not LR
class seq2fitness_handler:
    def __init__(self):
        
        model_path = 'nn-extrapolation-models/pretrained_models/other_models/gb1_lr_ckpt.pb'
        self.model = inference.restore_sess_from_pb(model_path)

        #print('model defined')

    def seq2fitness(self, seq):
        encoded_seq = encoding.encode(encoding="one_hot,aa_index", char_seqs=[seq])
        prediction = inference.run_inference(encoded_data=encoded_seq, sess=self.model)
        
        return prediction



