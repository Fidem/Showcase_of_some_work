import numpy as np


# This is a VERY trivial neural network example:
# Architecture is feed-forward 2 layer (input and output layer)
# Using gradient descent to optimise weights

def sigmoid(x):
    return 1 / (1 + np.exp(-x))


# Training input and output is arbitrary, but I want it to recognise that if the first
# input x1 = 1, then the output is more likely to be 1.

training_input = np.array([[1, 1, 0], [0, 0, 1], [0, 1, 0], [1, 1, 0], [1, 1, 1]])
training_output = np.array([[1, 0, 0, 1, 1]])
training_output = training_output.reshape(5, 1)  # Formatting output

np.random.seed(42)
# Randomly assigning initial values to weights and bias for symmetry breaking
# purposes.
weights = np.random.rand(3, 1)
bias = np.random.rand(1)
alpha = 0.05  # The Learning rate -> used in gradient descent calculation


# Derivative of the sigmoid function, also used in gradient descent/backpropogation
def sigmoid_der(x):
    return sigmoid(x) * (1 - sigmoid(x))


# Now we optimise the weights and bias!!
for epoch in range(20000):

    # Start with a feed forward to get values for each activation function, in this case
    # just the output
    # feedforward step1
    XW = np.dot(training_input, weights) + bias

    # feedforward step2
    z = sigmoid(XW)
    # Z here is the output.

    # Backpropagation step1
    # Calculate error (trivial for this simple neural network)
    error = z - training_output
    # Print the error within epochs so that we can view the error decreasing over time (hopefully)
    print(error.sum())

    # Backpropagation step 2
    dcost_dpred = error  # This is the difference between output and training value
    dpred_dz = sigmoid_der(z)  # Calculate the differential sigmoid(z) for backpropagation

    # This is the backpropagation error in the final layer.
    z_delta = dcost_dpred * dpred_dz

    inputs = training_input.T

    # Next are the iterative part of gradient descent.
    # We are moving the weights based on minimising the loss function (given by the
    # np.dot part of this equation)
    weights -= alpha * np.dot(inputs, z_delta)

    for num in z_delta:
        bias -= alpha * num

# The output shows the MSE loss function error and how it decreases over epochs


# This next bit just uses another arbitrary vector to see what the neural network will say
# about it.

single_point = np.array([1, 0, 0])
result = sigmoid(np.dot(single_point, weights) + bias)
print("Output Value: {}".format(result))

# This is perfect! the neural network is saying that the end value is VERY (0.99) likely
# to be 1.

# What about for the converse.
single_point = np.array([0, 1, 1])
result = sigmoid(np.dot(single_point, weights) + bias)
print("Output Value: {}".format(result))
# An absolutely tiny chance of being 1.  This is what I was trying to do!
# The neural network actually works!!!  On to bigger and better things now...
