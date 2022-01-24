import numpy as np
from ase.io import read, write
import pickle
import rascal
from rascal.neighbourlist.structure_manager import mask_center_atoms_by_species
from rascal.models import KRR

species  = [1,6,7,8] # list of chemical species for which to compute shieldings
compound = 'glycine' # compound, required for loading the correct ML shielding model
nmodels  = 16 # all committees are composed of 16 individual models
inname   = '/home/lumiaro/Documents/master_project/average_trajectories/glycine_test/PI_NVT/glycine_alpha/lmp_run/simulation.pos_0_ext.xyz' # filename/path of the input trajectory to calculate shieldings for
outname  = '/home/lumiaro/Documents/master_project/average_trajectories/glycine_test/PI_NVT/glycine_shieldings/pi_nvt_lmp_glycine_alpha_300K_00_w_cs1.xyz' # filename/path of the file to which the trajactory with shieldings should be written
outname_averaged  = '/home/lumiaro/Documents/master_project/average_trajectories/glycine_test/PI_NVT/glycine_shieldings/pi_nvt_lmp_glycine_alpha_300K_00_w_cs_averaged1.xyz'
nbatch   = 10000 # number of environments to be considered at a given time (to limit memory requirements)

def predict(f_test, soap, kern, feat, weights, alpha):

    x_test = soap.transform(f_test)
    print('features evaluated')

    y_pred   = KRR(weights, kern, feat, self_contributions=None, units = {'energy': 'ppm', 'length': 'AA'}).predict(x_test).reshape((-1, len(weights))).T
    print('models evaluated')
      
    # rescaling of differences and corresponding uncertainty estimates according to likelihood maximisation for obtained observed validation errors
    y_final  = y_pred * 0.0 + np.mean(y_pred, axis=0)
    y_final += alpha * (y_pred - np.mean(y_pred, axis=0))
    return y_final

def ave_batch_soap(f_test, soap):

    no_frames = len(f_test)

    x_test = soap.transform(f_test)
    print('features evaluated')

    x_test_array = x_test.get_features(soap)

    no_atoms_non_masked = np.sum(f_test[0].get_array("center_atoms_mask"))


    all_x_test_soaps = x_test_array.reshape(int(x_test_array.shape[0]/no_atoms_non_masked), no_atoms_non_masked , x_test_array.shape[1])

    sum_soap = np.sum(all_x_test_soaps, 0)

    return sum_soap, no_frames

def predict_averaged(kern, feat, weights, alpha, sum_soap, no_frames):

    soap_averaged = sum_soap/no_frames


    model = KRR(weights, kern, feat, self_contributions=None, units = {'energy': 'ppm', 'length': 'AA'})
   
    # also predict averages
    zeta = model.kernel._kwargs["zeta"]
    
    X = model.X_train.get_features()
    Y = soap_averaged
    XY = X.dot(Y.T)

    KNM = XY.T**zeta

    y_pred_averaged = np.dot(KNM, model.weights).reshape((-1)).reshape((-1, len(weights))).T

    y_final_averaged  = y_pred_averaged  * 0.0 + np.mean(y_pred_averaged , axis=0)
    y_final_averaged  += alpha * (y_pred_averaged  - np.mean(y_pred_averaged , axis=0))

    return y_final_averaged
    
traj = read(inname, ':')
# prepare test frames
bohr2ang = 0.529177

for ifrm, frm in enumerate(traj):
    # set arrays in which the chemical shieldings, uncertainty estimates, and individual predictions from the members of the committee will be stored
    frm.set_array('CS', np.zeros(frm.get_global_number_of_atoms()))
    frm.set_array('CSerr', np.zeros(frm.get_global_number_of_atoms()))
    frm.set_array('CSens', np.zeros((frm.get_global_number_of_atoms(), nmodels)))

    # convert cell parameters and atomic positions to Angstrom
    #frm.positions *= bohr2ang
    #frm.cell      *= bohr2ang

    # wrap all atoms into the unit cell
    frm.wrap(eps=1e-11)
print('frames loaded and prepared')

output_atom = traj[0].copy()
output_atom.positions = np.zeros(output_atom.positions.shape)
    
    
# run over chemical species
for sp in species:
    print('considering chemical species', sp)

    # load ShiftML model
    with open('/home/lumiaro/Documents/master_project/average_trajectories/glycine_test/glycine_models/' + compound + '_' + str(sp) + '.pickle','rb') as f:
        [soap, kern, feat, weig, alpha] = pickle.load(f)

    print('loaded shielding model')
    # mask all atomic centers except for species at hand
    nsp = []

    for ifrm, frm in enumerate(traj):

        frm.set_array('center_atoms_mask', np.zeros(frm.get_global_number_of_atoms(), dtype=bool))

        mask_center_atoms_by_species(frm, species_select=[sp,])

        nsp.append(np.sum(frm.get_array('center_atoms_mask')))

    output_atom.set_array('center_atoms_mask', np.zeros(frm.get_global_number_of_atoms(), dtype=bool))

    mask_center_atoms_by_species(output_atom, species_select=[sp,])


    # split trajectory into batches
    ibatch = [0]
    jbatch = []
    
    #Snbatch = np.min(nbatch, len(traj))

    for ib in range( int(np.sum(nsp)/nbatch + 0.5)):
        ifrm = 0
        msp = 0
        while msp < (ib + 1) * nbatch:
            msp = np.sum(nsp[:ifrm])
            ifrm += 1
            if ifrm == len(traj):
                ifrm += 1
                break
        jbatch.append(ifrm)
        ibatch.append(ifrm)
    print('split frames into', len(jbatch), 'batches')
    
    first_batch = True

    # predict for test frames
    for ib in range(int(np.sum(nsp)/nbatch + 0.5)):
        if first_batch:
            sum_soap, no_frames = ave_batch_soap(traj[ ibatch[ib] : jbatch[ib] ], soap)
            first_batch = False
        else:
            sum_soap_tmp, no_frames_tmp = ave_batch_soap(traj[ ibatch[ib] : jbatch[ib] ], soap)
            sum_soap += sum_soap_tmp
            no_frames += no_frames_tmp

        pred = predict(traj[ ibatch[ib] : jbatch[ib] ], soap, kern, feat, weig, alpha)

        print('predicted for batch', ib + 1)

        pred_mean = pred.mean(axis=0)
        pred_err  = pred.std(axis=0)


        # save trajectory with predicted CS for future use
        counter = {}
        counter[sp] = 0

        for ifrm, frm in enumerate( traj[ ibatch[ib] : jbatch[ib] ] ):

            # get CS
            cs = frm.get_array('CS')
            cserr = frm.get_array('CSerr')
            csens = frm.get_array('CSens')

            # set predicted CS for H and C
            mask = np.where(frm.get_atomic_numbers() == sp)[0]

            cs[mask] = pred_mean[counter[sp] : counter[sp] + len(mask)]
            cserr[mask] = pred_err[counter[sp] : counter[sp] + len(mask)]
            csens[mask] = pred.T[counter[sp] : counter[sp] + len(mask)]
            counter[sp] += len(mask)

            # set collected CS for frame
            frm.set_array('CS', cs)
            frm.set_array('CSerr', cserr)
            frm.set_array('CSens', csens)
            
        print('predictions stored for batch', ib+1)

    

    pred_averaged = predict_averaged(kern, feat, weig, alpha, sum_soap, no_frames)

    pred_mean_averaged = pred_averaged.mean(axis=0)
    pred_err_averaged  = pred_averaged.std(axis=0)

    cs = output_atom.get_array('CS')
    cserr = output_atom.get_array('CSerr')
    csens = output_atom.get_array('CSens')

    mask = np.where(frm.get_atomic_numbers() == sp)[0]

    cs[mask] = pred_mean_averaged[0 :  len(mask)]
    cserr[mask] = pred_err_averaged[0 :  len(mask)]
    csens[mask] = pred_averaged.T[0 :  len(mask)]

    # set collected CS for frame
    output_atom.set_array('CS', cs)
    output_atom.set_array('CSerr', cserr)
    output_atom.set_array('CSens', csens)

    write(outname_averaged, output_atom)
    print(outname_averaged)
    #break
    write(outname,traj)

