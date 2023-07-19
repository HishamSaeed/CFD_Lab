#include <string>
#include <was/blob.h>
#include <was/storage_account.h>

void upload_file(std::string dirName);

void delete_dir(std::string dirName, azure::storage::list_blob_item_iterator segment);