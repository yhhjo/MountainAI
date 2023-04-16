import torch 
import torch.nn as nn 
import math

"""Retrieved from PyTorch Documentation: https://pytorch.org/tutorials/beginner/transformer_tutorial.html"""

class PositionalEncoding(nn.Module):

    def __init__(self, d_model: int, dropout: float = 0.1, seqlen: int = 10):
        super().__init__()
        self.dropout = nn.Dropout(p=dropout)

        position = torch.arange(seqlen).unsqueeze(1)
        div_term = torch.exp(torch.arange(0, d_model, 2) * (-math.log(10000.0) / d_model))
        pe = torch.zeros(1, seqlen, d_model)
        pe[0, :, 0::2] = torch.sin(position * div_term)
        pe[0, :, 1::2] = torch.cos(position * div_term)
        self.register_buffer('pe', pe)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """
        Arguments:
            x: Tensor, shape ``[batch_size, seqlen, d_model]``
        """
        x = x + self.pe[:x.size(1)]
        return self.dropout(x)