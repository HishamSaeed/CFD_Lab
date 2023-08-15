#include <filesystem>
#include <cpprest/filestream.h>
#include <cpprest/containerstream.h>
#include "azure-blob-storage.hpp"

void upload_file(std::string dirName) {
    // Initialize storage account
    azure::storage::cloud_storage_account storage_account = azure::storage::cloud_storage_account::parse(parse_access_token());

    // Create a blob container
    azure::storage::cloud_blob_client blob_client = storage_account.create_cloud_blob_client();
    azure::storage::cloud_blob_container container = blob_client.get_container_reference(_XPLATSTR("vtpfilestoragecontainer"));
    
    // Return value is true if the container did not exist and was successfully created.
    container.create_if_not_exists();
    
    // Make the blob container publicly accessible
    azure::storage::blob_container_permissions permissions;
    permissions.set_public_access(azure::storage::blob_container_public_access_type::blob);
    container.upload_permissions(permissions);
    
    // List all blobs in the container

    azure::storage::list_blob_item_iterator segment = container.list_blobs();
    delete_dir(dirName, segment);

    // Upload a blob from a file
    for (const auto& entry : std::filesystem::directory_iterator("/home/hisham/dev/repos_test/CFD_Lab/Temp_Results/")) {
        concurrency::streams::istream input_stream = concurrency::streams::file_stream<uint8_t>::open_istream(_XPLATSTR(entry.path())).get();
        azure::storage::cloud_block_blob blob1 = container.get_block_blob_reference(_XPLATSTR(dirName + "/" + entry.path().filename().string()));
        blob1.upload_from_stream(input_stream);
        input_stream.close().wait();
    }
}

void delete_dir(std::string dirName, azure::storage::list_blob_item_iterator segment) {
    for (const azure::storage::list_blob_item& blobItem : segment) {
        if(!blobItem.is_blob()) {
            const azure::storage::list_blob_item_iterator dirSegment = blobItem.as_directory().list_blobs();
            delete_dir(dirName, dirSegment);
        } else {
            azure::storage::cloud_blob blob = blobItem.as_blob();
            if (blob.name().compare(0, dirName.length(), dirName) == 0) {
                ((azure::storage::cloud_blob)blob).delete_blob();
            }
        }
    }
}

std::string parse_access_token() {
    // Read the JSON file
    std::ifstream file("../credentials/access-token.json");
    if (!file.is_open()) {
        std::cerr << "Error opening JSON file." << std::endl;
        return NULL;
    }

    // Parse JSON data
    nlohmann::json json_data;
    file >> json_data;
    file.close();

    // Get the access token
    std::string access_token = json_data["access_token"];

    return access_token;
}