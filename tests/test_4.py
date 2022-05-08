import numpy as np
import SimpleFEM as simFem
import matplotlib.pyplot as plt
plt.interactive(True)

#Compare Simple FEM and results from https://scicomp.stackexchange.com/questions/35367/i-wrote-a-2d-finite-element-program-for-axial-loaded-plates-but-the-results-are

def test_4():

    # constants
    P = 18000  # Load (Newtons)
    W = 48  # Width of Beam (mm)
    H = 12  # Height of Beam (mm)
    Z = 1  # Thickness of Beam (mm)
    E_beam = 30E6  # Beam Elastic Modulus
    pr_beam = 0.45  # Poissons Ratio of the beam
    nds_x = 9  # number of nodes extending in the horizontal direction
    nds_y = 3  # number of nodes extending in the vertical direction

    nnds = nds_x * nds_y  # total number of nodes
    ndof = nnds * 2  # total number of degrees of freedom in the whole system
    nele = 2 * (nds_x - 1) * (nds_y - 1)  # total number of elements
    eper = 2 * (nds_x - 1)  # elements per element row

    ndcoor = np.zeros((nnds, 2))  # Table which stores the INITIAL coordinates (in terms of mm) for each node
    nd_rnc = np.zeros((nnds, 2))  # Table which stores the 'row and column' coordinates for each node
    nds_in_ele = np.zeros((nele, 3), dtype=np.int)  # the nodes which comprise each element
    glbStiff = np.zeros((ndof, ndof))  # global stiffness matrix (GSM)
    lst_wallnds = []  # List of nodes (indices) which are coincident with the rigid wall on the left
    lst_wallnds.clear()
    lst_walldofs = []  # dofs indices of nodes coincident with the rigid wall
    lst_walldofs.clear()
    lst_endnds = []  # nodes on the free edge of the beam
    lst_endnds.clear()

    nnf_msg = 'Node not found!'


    # Function 'node_by_rnc' returns the index of the node which has the same row and column as the ones input (in_row, in_col)
    def node_by_rnc(in_row, in_col, start_mrow):  # 'start_mrow' == where the func starts searching (for efficiency)
        run = True
        row = start_mrow
        while run == True:
            if row > nnds - 1:
                run = False
            elif nd_rnc[row][0] == in_row and nd_rnc[row][1] == in_col:
                run = False
            else:
                row = row + 1
        if row > nnds - 1:
            return nnf_msg  # returns error message
        else:
            return row


    # Function 'add_to_glbStiff' takes a local stiffness matrix and adds the value of each 'cell' to the corrosponding cell in the GSM
    def add_to_glbStiff(in_mtrx, nd_a, nd_b, nd_c):
        global glbStiff
        # First column in local stiffness matrix (LSM) is the x-DOF of Node A, second is the y-DOF of Node A, third is the x-DOF of Node B, etc. (same system for rows; the matrix is symmetric)
        dofs = [2 * nd_a, 2 * nd_a + 1, 2 * nd_b, 2 * nd_b + 1, 2 * nd_c,
                2 * nd_c + 1]  # x-DOF for a node == 2*[index of the node]; y-DOF for node == 2*[node index]+1
        for r in range(0, 6):  # LSMs are always 6x6
            for c in range(0, 6):
                gr = dofs[r]  # gr == row in global stiffness matrix
                gc = dofs[c]  # gc == column in global stiffness matrix
                glbStiff[gr][gc] = glbStiff[gr][gc] + in_mtrx[r][
                    c]  # Add the value of the LSM 'cell' to what's already in the corrosponding GSM cell


    for n in range(0, nnds):  # puts node coordinates and rnc indices into matrix
        row = n // nds_x
        col = n % nds_x
        nd_rnc[n][0] = row
        nd_rnc[n][1] = col
        ndcoor[n][0] = col * (W / (nds_x - 1))
        ndcoor[n][1] = row * (H / (nds_y - 1))

        if col == 0:
            lst_wallnds.append(n)
        elif col == nds_x - 1:
            lst_endnds.append(n)

    for e in range(0, nele):  # FOR EVERY ELEMENT IN THE SYSTEM...
        # ...DETERMINE NODES WHICH COMPRISE THE ELEMENT
        erow = e // eper  # erow == the row which element 'e' is on
        eor = e % eper  # element number on row (i.e. eor==0 means the element is attached to rigid wall)
        if eor % 2 == 0:  # downwards-facing triangle
            nd_a_col = eor / 2
            nd_b_col = eor / 2
            nd_c_col = (eor / 2) + 1
            nd_a = node_by_rnc(erow, nd_a_col, nds_x * erow)
            nd_b = node_by_rnc(erow + 1, nd_b_col, nds_x * erow)
            nd_c = node_by_rnc(erow, nd_c_col, nds_x * erow)
        else:  # upwards-facing triangle
            nd_a_col = (eor // 2) + 1
            nd_b_col = (eor // 2) + 1
            nd_c_col = eor // 2
            nd_a = node_by_rnc(erow + 1, nd_a_col, nds_x * (erow + 1))
            nd_b = node_by_rnc(erow, nd_b_col, nds_x * erow)
            nd_c = node_by_rnc(erow + 1, nd_c_col, nds_x * (erow + 1))
        if nd_a != nnf_msg and nd_b != nnf_msg and nd_c != nnf_msg:  # assign matrix element values if no error
            nds_in_ele[e][0] = int(nd_a)
            nds_in_ele[e][1] = int(nd_b)
            nds_in_ele[e][2] = int(nd_c)
        else:  # raise error
            print(nnf_msg)

    '''
       MY CODE
    '''
    nodes = np.array([ ndcoor[:, 0], ndcoor[:, 1] ])
    connectivity = nds_in_ele



    alumininium = simFem.Elastic(30E6, 0.45)
    t = 1.0
    boundaryConditions = {
        "fixed": {0: "all", 9: "all", 18: "all"},
        #"applied force": {4: [0, 6.67 * 1000], 9: [0, 6.67 * 1000], 26: [0, 26.67 * 1000]}
        "applied force": {8: [6000, 0], 17: [6000, 0], 26: [6000, 0]}

    }
    plate = simFem.Body(simFem.createTriangles(connectivity, nodes, alumininium, t), np.shape(nodes)[1])
    plate.assembleGlobalStiffness()
    with open('MyStiff.npy', 'wb') as f:
        np.save(f, plate.K)
    plate.applyBcs(boundaryConditions)
    plate.solve()

    fig, ax = plt.subplots()
    ax.scatter(nodes[0, :], nodes[1, :])
    for i in range(np.shape(nodes)[1]):
        ax.annotate(i, (nodes[0, i], nodes[1, i]))
    plt.triplot(nodes[0, :], nodes[1, :], connectivity)
    plt.show()



