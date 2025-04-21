from stellargraph.mapper import  CorruptedGenerator, HinSAGENodeGenerator
from sklearn.preprocessing import StandardScaler
from stellargraph import StellarGraph
from stellargraph.utils import plot_history
from stellargraph.layer import DeepGraphInfomax, HinSAGE

import pandas as pd

import tensorflow as tf
from keras.api._v2.keras import Model
from keras.api._v2.keras.optimizers import Adam
from keras.api._v2.keras.callbacks import EarlyStopping

tf.random.set_seed(42)

epochs = 1000
# Early stopping to prevent overfitting
es = EarlyStopping(monitor="loss", min_delta=0.001, patience=10, verbose=1, mode='min')

# Function to run Deep Graph Infomax
def run_deep_graph_infomax(
    base_model, generator, epochs, reorder=lambda sequence, subjects: subjects
):
    lncRNA_nodes = LPPI_graph.nodes(node_type='lncRNA')

    corrupted_generator = CorruptedGenerator(generator)
    lncRNA_gen = corrupted_generator.flow(lncRNA_nodes)
    infomax = DeepGraphInfomax(base_model, corrupted_generator)

    x_in, x_out = infomax.in_out_tensors()

    model = Model(inputs=x_in, outputs=x_out)
    model.compile(loss=tf.nn.sigmoid_cross_entropy_with_logits, optimizer=Adam(learning_rate=1e-2))
    history = model.fit(lncRNA_gen, epochs=epochs, verbose=2, callbacks=[es])
    plot_history(history)
    x_emb_in, x_emb_out = base_model.in_out_tensors()

    emb_model = Model(inputs=x_emb_in, outputs=x_emb_out)
    node_gen = generator.flow(lncRNA_nodes, targets=None)
    embedding = emb_model.predict(node_gen)

    return embedding

# Tuning hyperparameters

#proteins = pd.read_csv(f"../../Annotate/mouse/transformed_protein_annotation.csv")
#lncRNAs = pd.read_csv(f"../../Annotate/mouse/valid_heart_annotation.csv")
#LPPI = pd.read_csv(f'../../Annotate/mouse/weighted_valid_inter.csv')

#proteins.set_index(proteins['protein_ID'], inplace=True)
#proteins = proteins.drop('protein_ID', axis=1)
#lncRNAs.set_index(lncRNAs['lncRNA_ID'], inplace=True)
#lncRNAs = lncRNAs.drop('lncRNA_ID', axis=1)

#scaler_lncRNAs = StandardScaler()
#lncRNAs_scaled = scaler_lncRNAs.fit_transform(lncRNAs)

## create StellarGraph
#LPPI_graph = StellarGraph({"lncRNA": pd.DataFrame(lncRNAs_scaled, index=lncRNAs.index),
#                        "protein": pd.DataFrame(proteins, index=proteins.index)},
#                        LPPI)
#print(LPPI_graph.info())

#layer_size_1=[16,32,64]
#layer_size_2=[32,64,128,256]

#samples_num_1=[5,10,15,20]
#samples_num_2=[10,15,20,25]

#for i in samples_num_1:
#    for j in samples_num_2:
#        hinsage_generator = HinSAGENodeGenerator(
#            LPPI_graph, batch_size=32, num_samples=[i,j], head_node_type="lncRNA"
#        )
#        hinsage_model = HinSAGE(
#            layer_sizes=[32,64], activations=["LeakyReLU", "linear"], generator=hinsage_generator)
#        print(f"Train begin: {i} for first layer,{j} for second layer.")
#        node_embedding = run_deep_graph_infomax(hinsage_model, hinsage_generator, epochs=epochs)

#        embedding_df = pd.DataFrame(node_embedding, index=LPPI_graph.nodes(node_type='lncRNA'))
#        embedding_df.to_csv(f'mouse/samples_num/lncRNA_embeddings_MEL_{i}_{j}.csv', header=None)

human_tissues=['heart','lung', 'stomach']
mouse_tissues=['heart','lung', 'brain']

# Modify this parameter to train models for different species.
species = 'mouse'

# Iterate over the list of tissues corresponding to the given species.
for tissue in mouse_tissues:

    proteins = pd.read_csv(f"../../Annotate/{species}/transformed_protein_annotation.csv")
    lncRNAs = pd.read_csv(f"../../Annotate/{species}/valid_{tissue}_annotation.csv")
    LPPI = pd.read_csv(f'../../Annotate/{species}/weighted_valid_inter.csv')

    proteins.set_index(proteins['protein_ID'], inplace=True)
    proteins = proteins.drop('protein_ID', axis=1)
    lncRNAs.set_index(lncRNAs['lncRNA_ID'], inplace=True)
    lncRNAs = lncRNAs.drop('lncRNA_ID', axis=1)

    # Normalize the lncRNA features
    scaler_lncRNAs = StandardScaler()
    lncRNAs_scaled = scaler_lncRNAs.fit_transform(lncRNAs)

    # create StellarGraph
    LPPI_graph = StellarGraph({"lncRNA": pd.DataFrame(lncRNAs_scaled, index=lncRNAs.index),
                            "protein": pd.DataFrame(proteins, index=proteins.index)},
                            LPPI)
    print(LPPI_graph.info())

    hinsage_generator = HinSAGENodeGenerator(
        LPPI_graph, batch_size=32, num_samples=[20,20], head_node_type="lncRNA"
    )
    hinsage_model = HinSAGE(
        layer_sizes=[16,32], activations=["LeakyReLU", "linear"], generator=hinsage_generator)
    
    node_embedding = run_deep_graph_infomax(hinsage_model, hinsage_generator, epochs=epochs)

    # Save the embeddings to a CSV file
    embedding_df = pd.DataFrame(node_embedding, index=LPPI_graph.nodes(node_type='lncRNA'))
    embedding_df.to_csv(f'./{species}/lncRNA_embeddings_{tissue}.csv', header=None)