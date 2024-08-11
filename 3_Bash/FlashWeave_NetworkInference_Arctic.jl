# Load FlashWeave
using FlashWeave

# Define input path for the count data
count_data_path = string(""/home/alle/WORKING_DIR/Cristina/Networks/FlashWeaveNetwork_Poles_CountData.csv")

# Define input path for the metadata
meta_data_path = string("home/alle/WORKING_DIR/Cristina/Networks/FlashWeaveNetwork_Poles_Metadata.csv")

# Calculate network
network = learn_network(count_data_path, meta_data_path, sensitive=true, heterogeneous=false)

# Save network
save_network(string("/home/alle/WORKING_DIR/Cristina/Networks/NetworkFlashWeave.edgelist"), network)
