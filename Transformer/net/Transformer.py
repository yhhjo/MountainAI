import torch
import torch.nn as nn
from Encoder import EncoderBlock
from Decoder import DecoderBlock

class Transformer(nn.Module):

    def __init__(self, d_input: int, seqlen: int, d_output: int, d_model: int, q: int, v: int, h: int, Ne: int, Nd: int, batch_size: int, dropout = 0.1):
        """Transformer Model from Attention is All You Need, d2l.ai, and Cohen et. all 2022. 
        The transformer consists of an embedding layer, a stack of N encoding blocks, a stack of N decoding blocks,
        followed by a linear transformation and an activation function, 

        Parameters 
        ----------
        d_input: 
            Dimension of input vector
        d_output:
            Dimension output
        d_model:
            Dimension model's latent space
        q:
            Query dimensionality
        v:
            Value dimensionality
        h: 
            Number of heads 
        Ne:
            Number of encoding layers
        Nd: 
            Number of decoding layers
        dropout:
            Dropout rate for dot product attention module.
        """

        super().__init__()
        
        self._batch_size = batch_size
        self._encodingStack = []
        self._decodingStack = []
        self._embeddingLayer = nn.Linear(d_input, d_model)
        self._out1 = nn.Linear(d_model, d_output)

        for i in range(Ne):
            self._encodingStack.append(EncoderBlock(d_model, q, v, h, dropout))
        for i in range(Nd):
            self._decodingStack.append(DecoderBlock(d_model, q, v, h, dropout))

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """
        Parameters: 
        ----------
        x:
            Tensor of shape (batch_size, seqlen, d_output)

        Returns
        ----------
        x:
            Tensor shape (batch_size, d_output)
        """
        x = self._embeddingLayer(x)

        for e in self._encodingStack:
            x = e(x)

        
        context = x

        for d in self._decodingStack:
            x = d(x, context)

        x = x.mean(dim=1)  # Global mean pooling along the `seqlen` dimension
        x = self._out1(x)
        return x