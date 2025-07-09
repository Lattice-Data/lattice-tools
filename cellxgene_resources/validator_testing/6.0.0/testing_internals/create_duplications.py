import anndata as ad
import numpy as np
import dask.array as da
import scipy.sparse as sp


def create_duplications(adata:ad.AnnData,barcodes:list[str],n=int):
        '''
        Duplicates raw cell counts within a valid h5ad by choosing random rows and replacing them with a specified index.
        barcodes = list of indices to duplicate; n = number of times to duplicate each barcode

        :rtype Raw Anndata with duplications.
        '''
        list_for_obs = []
        total_needed = n * len(barcodes)
        if total_needed > len(adata.obs):
                raise ValueError("Not enough rows in the DataFrame to sample without replacement.")

        candidate_indices = adata.obs.index[~adata.obs.index.isin(barcodes)].tolist()
        if "in_tissue" in adata.obs.columns:
                obs_to_keep_in_tissue = set(adata.obs[adata.obs["in_tissue"]== 1].index)
                candidate_indices = list(set(candidate_indices) - obs_to_keep_in_tissue)

        random_indices = np.random.choice(candidate_indices, size=total_needed, replace=False)  # Sample total unique indices minus those specified in barcodes list
        random_indices_list = [random_indices[i * n:(i + 1) * n].tolist() for i in range(len(barcodes))]
        obs_name_to_index = {name: i for i, name in enumerate(adata.obs_names)}
        for barcode, target_barcodes in zip(barcodes, random_indices_list):
                obs_indx = obs_name_to_index[barcode]
                obs_meta = adata.obs.iloc[obs_indx].copy()
                list_for_obs.append([barcode,obs_indx,obs_meta,target_barcodes])

        dict_for_mtx = {}
        count = 0
        for barcode,obs_indx,obs_meta,target_barcodes in list_for_obs:
                for tb in target_barcodes:
                        adata.obs.loc[tb] = obs_meta
                target_indices = [obs_name_to_index[name] for name in target_barcodes]
                dict_for_mtx[f"indices_{count}"] = (obs_indx,target_indices)
                count += 1

        if adata.raw:
                dict_for_mtx["raw_adata"] = adata.raw
        else:
                dict_for_mtx["raw_adata"] = adata

        raw_adata = dict_for_mtx["raw_adata"]
        raw_matrix = dict_for_mtx["raw_adata"].X
        is_dask = hasattr(raw_matrix, "compute")
        if is_dask:
                raw_matrix = raw_matrix.compute().toarray()
                for key,value in dict_for_mtx.items():
                        if "indices" in key:
                                for i in value[1]:
                                        raw_matrix[i,:] = raw_matrix[value[0], :]

                raw_matrix_csr = sp.csr_matrix(raw_matrix)
                dup_raw_adata = ad.AnnData(raw_matrix_csr, obs=adata.obs, var=raw_adata.var, uns=adata.uns, obsm=adata.obsm)
                dup_raw_adata.X = da.from_array(dup_raw_adata.X,chunks=raw_adata.X.chunks)

        else:
                raw_matrix = raw_matrix.toarray()
                for key,value in dict_for_mtx.items():
                        if "indices" in key:
                                for i in value[1]:
                                        raw_matrix[i,:] = raw_matrix[value[0], :]

                dup_raw_adata = ad.AnnData(raw_matrix, obs=adata.obs, var=raw_adata.var, uns=adata.uns, obsm=adata.obsm)
        return dup_raw_adata
