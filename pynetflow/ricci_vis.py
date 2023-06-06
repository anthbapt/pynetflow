from pynetflow.ricci_curvature import forman_degree, forman_hyperedge
import matplotlib.pyplot as plt
from matplotlib import cm
import networkx as nx
import pandas as pd
import xgi
import cv2
import re
import os


def visualize_flow(flow, path = None) :
    if path == None:
        path = os.path.dirname(os.path.realpath(__file__))
        os.chdir(path)
    newfolder = 'ricci_flow'
    if not os.path.exists(newfolder) :
        os.makedirs(newfolder)
    os.chdir(path + newfolder)

    data = flow['df']
    G_list = [data[data[0]==i][[1,3]] for i in range(n_iter)]
    G_nodes = list()
    G_net = list()
    for k in range(n_iter) :
        G_temp = G_list[k]
        G_temp[[0,1]] = pd.DataFrame(G_temp[1].tolist(), index=G_temp.index)
        G_net.append(nx.from_pandas_edgelist(G_temp, 0, 1, [3]))
        temp = list()
        for i in range(len(G.nodes)) :
            temp.append(sum(G_list[k][G_list[k][0] == i][3]))
        G_nodes.append(temp)
        
    
    pos = nx.spring_layout(G_net[0])
    for k in range(n_iter) :
        edges,weights = zip(*nx.get_edge_attributes(G_net[k],3).items())
        nx.draw(G_net[k], pos, node_color='b', edgelist=edges, edge_color=weights, width=2.0, edge_cmap=plt.cm.seismic)
        plt.savefig('net_ricci_curvature_' + str(k) + '.png', format = 'png', dpi = 500)
        plt.close()
        
        
    for k in range(n_iter) :
        nx.draw(G_net[k], pos, nodelist = sorted(list(G_net[k].nodes)), node_color=G_nodes[k], cmap=plt.cm.seismic)
        plt.savefig('net_ricci_curvature_nodes' + str(k) + '.png', format = 'png', dpi = 500)
        plt.close()
        
    return None


def make_video(path = None, image_folder = 'ricci_flow', video_name = 'ricci_flow.avi', fps = 10):
    if path == None:
        path = os.path.dirname(os.path.realpath(__file__))
        os.chdir(path)
    os.chdir(path)
    images = [img for img in os.listdir(image_folder) if img.endswith(".png")]
    images = sorted(images, key=lambda s: int(re.search(r'\d+', s).group()))
    frame = cv2.imread(os.path.join(image_folder, images[0]))
    height, width, layers = frame.shape
    resolution = (width,height)
    codec = cv2.VideoWriter_fourcc(*'x264')
    video = cv2.VideoWriter(video_name, codec, fps, resolution)
    
    for image in images:
        video.write(cv2.imread(os.path.join(image_folder, image)))
    
    cv2.destroyAllWindows()
    video.release()
    
    return None


def visualize_order(SC, p, path = None, save = True):
    if path == None:
        path = os.path.dirname(os.path.realpath(__file__))
        os.chdir(path)
    pos = xgi.barycenter_spring_layout(SC)
    if p == 0 :
        vmin = min(forman_degree(SC, SC.nodes).values())
        vmax = max(forman_degree(SC, SC.nodes).values())
        plt.figure(figsize=(10, 10))
        ax = plt.subplot(111)
        xgi.draw(
            SC,
            pos,
            node_fc=SC.nodes.forman_degree,
            node_fc_cmap=cm.seismic, edge_fc = 'white',
            ax=ax)
        sm = plt.cm.ScalarMappable(cmap=cm.seismic, norm=plt.Normalize(vmin = vmin, vmax=vmax))
        sm._A = []
        plt.colorbar(sm)
        if save == True:
            plt.savefig('test_ricci_curvature_nodes.png', format = 'png', dpi = 500)
            plt.close()
    elif p == 1 :
        SC_sub = xgi.SimplicialComplex()
        SC_sub.add_simplices_from(SC.edges.filterby('order', p).members())
        SC_sub.add_nodes_from(SC.nodes - SC_sub.nodes)
        vmin = min(forman_hyperedge(SC, SC.nodes, p).values())
        vmax = max(forman_hyperedge(SC, SC.nodes, p).values())
        plt.figure(figsize=(10, 10))
        ax = plt.subplot(111)
        xgi.draw(
            SC_sub,
            pos,
            dyad_color = SC.edges.filterby('order', p).forman_hyperedge(p),
            dyad_color_cmap = cm.seismic, edge_fc = 'white',
            ax=ax)
        sm = plt.cm.ScalarMappable(cmap=cm.seismic, norm=plt.Normalize(vmin = vmin, vmax=vmax))
        sm._A = []
        plt.colorbar(sm)
        if save == True:
            plt.savefig('test_ricci_curvature_edges.png', format = 'png', dpi = 500)
            plt.close()
    elif p == 2 :
        # triangles
        SC_sub = xgi.SimplicialComplex()
        SC_sub.add_simplices_from(SC.edges.filterby('order', p).members())
        SC_sub.add_nodes_from(SC.nodes - SC_sub.nodes)
        vmin = min(forman_hyperedge(SC, SC.nodes, p).values())
        vmax = max(forman_hyperedge(SC, SC.nodes, p).values())
        plt.figure(figsize=(10, 10))
        ax = plt.subplot(111)
        xgi.draw(
            SC_sub,
            pos,
            edge_fc = SC.edges.filterby('order',p).forman_hyperedge(p),
            edge_fc_cmap = cm.seismic,
            ax=ax)
        sm = plt.cm.ScalarMappable(cmap=cm.seismic, norm=plt.Normalize(vmin = vmin, vmax=vmax))
        sm._A = []
        plt.colorbar(sm)
        if save == True:
            plt.savefig('test_ricci_curvature_triangles.png', format = 'png', dpi = 500)
            plt.close()
    else :
        SC_sub = xgi.SimplicialComplex()
        SC_sub.add_simplices_from(SC.edges.filterby('order', p).members())
        SC_sub.add_nodes_from(SC.nodes - SC_sub.nodes)
        vmin = min(forman_hyperedge(SC, SC.nodes, p).values())
        vmax = max(forman_hyperedge(SC, SC.nodes, p).values())
        plt.figure(figsize=(10, 10))
        ax = plt.subplot(111)
        xgi.draw(
            SC_sub,
            pos,
            edge_fc = SC.edges.filterby('order',p).forman_hyperedge(p),
            edge_fc_cmap = cm.seismic,
            ax=ax)
        sm = plt.cm.ScalarMappable(cmap=cm.seismic, norm=plt.Normalize(vmin = vmin, vmax=vmax))
        sm._A = []
        plt.colorbar(sm)
        if save == True:
            plt.savefig('test_ricci_curvature_hyperedges' + str(p) + '.png', format = 'png', dpi = 500)
            plt.close()
        
    return None

    
def visualize_all(SC, path = None, save = True):
    if path == None:
        path = os.path.dirname(os.path.realpath(__file__))
        os.chdir(path)
    pos = xgi.barycenter_spring_layout(SC)
    max_order = len(max(SC.edges.members(), key = len))
    # max_order-1 to avoid the issue detailed below
    for k in range(max_order) :
        if k == 0 :
            vmin = min(forman_degree(SC, SC.nodes).values())
            vmax = max(forman_degree(SC, SC.nodes).values())
            plt.figure(figsize=(10, 10))
            ax = plt.subplot(111)
            xgi.draw(
                SC,
                pos,
                node_fc=SC.nodes.forman_degree,
                node_fc_cmap=cm.seismic, edge_fc = 'white',
                ax=ax)
            sm = plt.cm.ScalarMappable(cmap=cm.seismic, norm=plt.Normalize(vmin = vmin, vmax=vmax))
            sm._A = []
            plt.colorbar(sm)
            if save == True:
                plt.savefig('test_ricci_curvature_nodes.png', format = 'png', dpi = 500)
                plt.close()
        elif k == 1 :
            SC_sub = xgi.SimplicialComplex()
            SC_sub.add_simplices_from(SC.edges.filterby('order', k).members())
            SC_sub.add_nodes_from(SC.nodes - SC_sub.nodes)
            vmin = min(forman_hyperedge(SC, SC.nodes, k).values())
            vmax = max(forman_hyperedge(SC, SC.nodes, k).values())
            plt.figure(figsize=(10, 10))
            ax = plt.subplot(111)
            xgi.draw(
                SC_sub,
                pos,
                dyad_color = SC.edges.filterby('order', k).forman_hyperedge(k),
                dyad_color_cmap = cm.seismic, edge_fc = 'white',
                ax=ax)
            sm = plt.cm.ScalarMappable(cmap=cm.seismic, norm=plt.Normalize(vmin = vmin, vmax=vmax))
            sm._A = []
            plt.colorbar(sm)
            if save == True:
                plt.savefig('test_ricci_curvature_edges.png', format = 'png', dpi = 500)
                plt.close()
        elif k == 2 :
            # triangles
            SC_sub = xgi.SimplicialComplex()
            SC_sub.add_simplices_from(SC.edges.filterby('order', k).members())
            SC_sub.add_nodes_from(SC.nodes - SC_sub.nodes)
            vmin = min(forman_hyperedge(SC, SC.nodes, k).values())
            vmax = max(forman_hyperedge(SC, SC.nodes, k).values())
            plt.figure(figsize=(10, 10))
            ax = plt.subplot(111)
            xgi.draw(
                SC_sub,
                pos,
                edge_fc = SC.edges.filterby('order',k).forman_hyperedge(k),
                edge_fc_cmap = cm.seismic,
                ax=ax)
            sm = plt.cm.ScalarMappable(cmap=cm.seismic, norm=plt.Normalize(vmin = vmin, vmax=vmax))
            sm._A = []
            plt.colorbar(sm)
            if save == True:
                plt.savefig('test_ricci_curvature_triangles.png', format = 'png', dpi = 500)
                plt.close()
        else :
            SC_sub = xgi.SimplicialComplex()
            SC_sub.add_simplices_from(SC.edges.filterby('order', k).members())
            SC_sub.add_nodes_from(SC.nodes - SC_sub.nodes)
            vmin = min(forman_hyperedge(SC, SC.nodes, k).values())
            vmax = max(forman_hyperedge(SC, SC.nodes, k).values())
            plt.figure(figsize=(10, 10))
            ax = plt.subplot(111)
            xgi.draw(
                SC_sub,
                pos,
                edge_fc = SC.edges.filterby('order',k).forman_hyperedge(k),
                edge_fc_cmap = cm.seismic,
                ax=ax)
            sm = plt.cm.ScalarMappable(cmap=cm.seismic, norm=plt.Normalize(vmin = vmin, vmax=vmax))
            sm._A = []
            plt.colorbar(sm)
            if save == True:
                plt.savefig('test_ricci_curvature_hyperedges' + str(k) + '.png', format = 'png', dpi = 500)
                plt.close()
            
    return None