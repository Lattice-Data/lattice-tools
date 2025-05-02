import anndata as ad
import random

def create_duplications(adata:ad.AnnData,i:str,n=int):
        '''
        Duplicates raw cell counts within a valid h5ad by choosing a number of random rows and replacing them with a specified index.
        i = index to duplicate; n = number of times to duplicate

        :rtype Raw Anndata with duplications.
        '''
        print(adata.shape)
        obs_loc = adata.obs.index.get_loc(i)
        obs_row_for_dup = adata.obs.iloc[obs_loc].copy()
        random_indices = [random.randint(1, adata.obs.shape[0]-1) for m in range(n)]
        print(f"Replacing rows:{random_indices} with row {obs_loc}")
        for o_indx in random_indices:
            print(o_indx)
            adata.obs.iloc[o_indx] = obs_row_for_dup.values

        to_validate=[]

        if adata.raw:
                print("Raw counts found in adata.raw")
                to_validate.append(adata.raw)
        else:
                print("Raw counts found in adata")
                to_validate.append(adata)

        for raw_adata in to_validate:
                raw_matrix = raw_adata.X
                raw_matrix_row = raw_matrix[obs_loc, :]
                for m_indx in random_indices:
                        raw_matrix[m_indx] = raw_matrix_row

                raw_adata = ad.AnnData(raw_matrix, obs=adata.obs, var=adata.var, uns=adata.uns, obsm=adata.obsm)
                return raw_adata