module common_variables_RF

    implicit none

    character(len=100) :: results_folder_name, log_folder_name
    character(len=20) , parameter :: results_path = "./results"
    character (len=30), parameter :: log_filename = "log_proc_"
    character(len=20) , parameter :: xStep_folder_name = "xStep"
    character(len=20) , parameter :: Nmc_folder_name = "Nmc"
    character(len=20) , parameter :: corrL_folder_name = "corrL"
    character(len=200) :: single_path
    character(len=200) :: BBoxFileName = "BBox"


end module common_variables_RF
