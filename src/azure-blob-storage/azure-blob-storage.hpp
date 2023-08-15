#include <string>
#include <iostream>
#include <fstream>
#include <was/blob.h>
#include <was/storage_account.h>
#include "../json.hpp"

void upload_file(std::string dirName);

void delete_dir(std::string dirName, azure::storage::list_blob_item_iterator segment);

std::string parse_access_token();