# Import necessary libraries here
from os import listdir
from os.path import isfile, join, exists
import numpy as np
from simple_inference import encoding
from simple_inference import inference_lr as inference  # needed to remove training_ph from inference

class seq2fitness_handler:
    def __init__(self):
        
        model_path = 'pretrained_models/other_models/gb1_lr_ckpt.pb'
        self.model = inference.restore_sess_from_pb(model_path)

        #print('model defined')

    def seq2fitness(self, seq):
        encoded_seq = encoding.encode(encoding="one_hot,aa_index", char_seqs=[seq])
        prediction = inference.run_inference(encoded_data=encoded_seq, sess=self.model)
        
        return prediction



